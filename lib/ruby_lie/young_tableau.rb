module RubyLie

  class YoungTableau
    include Enumerable

    # The element of each row is a specific node from a vector representation
    attr_reader :rows
    attr_reader :vector_rep

    def initialize(params)
      if params[:vector_rep]
        @vector_rep = params[:vector_rep]
      else
        @vector_rep = nil
      end

      if params[:rows]
        @rows = Array.new

        params[:rows].each_with_index do |row, index|
          @rows[index] = Array.new
          row.each_with_index do |node, index2|
            @rows[index][index2] = node
          end
        end
      elsif params[:partition]
        partition = params[:partition]
        @rows = Array.new

        partition.each_with_index do |num, index|
          # Create new array at this row with all nil for now
          @rows[index] = Array.new(num, nil)
        end
      end

    end

    def partition
      @rows.each_with_index.map do |row, index|
        row.size
      end
    end

    def self.from_highest_weight(highest_weight)
      # When spinor rep, don't use
      # TODO make it so when even, still use YoungTableau,
      # e.g. D_4, 2*omega_3 rep has column [h_1,h_2,h_3,h_4]
      case highest_weight.algebra.alg
      when :alg_B
        return nil if ((highest_weight * highest_weight.algebra.alpha_dual(highest_weight.algebra.rank)) % 2) == 1
      when :alg_D
        return nil if ((highest_weight * highest_weight.algebra.alpha_dual(highest_weight.algebra.rank-1)) % 2) == 1
        return nil if ((highest_weight * highest_weight.algebra.alpha_dual(highest_weight.algebra.rank)) % 2) == 1
      end

      # First get partition
      partition = (1..highest_weight.algebra.rank).inject([]) do |all, i|
        blocks_in_row = (i..highest_weight.algebra.rank).inject(0) do |sum, cur|
          sum + highest_weight.algebra.alpha_dual(cur)*highest_weight
        end
        # Only accept the upper non-zero size rows
        if blocks_in_row == 0
          # Return the final partition array
          break(all)
        end
        all << blocks_in_row
      end

      tableau = YoungTableau.new(:vector_rep => highest_weight.algebra.vector_rep, :partition => partition)
      tableau.init_with_vector_rep(highest_weight.algebra.vector_rep)

      return tableau
    end

    def init_with_vector_rep(vector_rep)
      @vector_rep = vector_rep
      @rows.each_with_index do |row, index|
        @rows[index] = row.map {|i| vector_rep.node_at_index(index)}
      end
    end

    def each
      self.each_with_index do |elem, i,j|
        yield elem
      end
    end

    def each_with_index
      @rows.each_with_index do |row, index|
        i = row.size-1
        loop do
          yield row[i], index, i
          i -= 1
          break if i < 0
        end
      end
    end
    
    def each_col_at(i)
      (0..@rows.size-1).each do |row|
        break if rows[row][i].nil?
        yield rows[row][i], row
      end
    end

    def each_col_at_backwards(i)
      row = @rows.size-1
      loop do
        yield rows[row][i], row if not rows[row][i].nil?
        break if row == 0
        row -= 1
      end
    end


    def ==(tableau)
      case tableau
      when YoungTableau
        return self.rows == tableau.rows
      when Array
        return self.rows == tableau
      else
        return false
        #throw TypeError, "#{tableau.class} must be of type YoungTableau or Array"
      end

      return false
    end

    def to_s
      s = ""
      @rows.each_with_index do |row, index|
        row.each_with_index do |node, index2|
          s += "#{@vector_rep.node_to_index_hash[node]} "
        end
        s += "\n"
      end
      return s
    end
        

    # Given hash of ways to increment each element of the Young tableau,
    # increment one such block going from right to left, top to bottom
    def next_tableau(root_index)
      if @vector_rep.nil? or root_index > @vector_rep.algebra.rank or root_index < 0
        return nil
      end
      
      if root_index > 0
        col = rows[0].size-1
        
        loop do
          tableau_test = try_col_with_root_index(col, root_index)
          return tableau_test if not tableau_test.nil?
          break if col <= 0
          col -= 1
        end
      elsif root_index == 0
        (0..rows[0].size-1).each do |col|
          tableau_test = try_col_with_root_index(col, root_index)
          return tableau_test if not tableau_test.nil?
        end
      end

      return nil
    end

    # Try to increment the current column with simple root "root_index"
    def try_col_with_root_index(col, root_index)

      # This is the procedure for checking a given node element
      # at the given row
      check_elem_of_row = Proc.new do |elem, row|
        @vector_rep.each_chain_link(root_index) do |from, to|
          if from == elem
            tableau = try_next_tableau(row, col, to, root_index)
            return tableau if not tableau.nil?
          end
        end
      end

      if root_index > 0
        self.each_col_at(col) do |elem, row|
          check_elem_of_row.call(elem,row)
        end

      elsif root_index == 0
        self.each_col_at_backwards(col) do |elem, row|
          check_elem_of_row.call(elem,row)
        end
      end

      return nil
    end


    # Try to increment the current column/row to vector "to"
    # from simple root "root_index"
    def try_next_tableau(row, col, to, root_index)
      if root_index > 0
        return nil if not verify_row_col_to(row, col, to)

        tableau = YoungTableau.new(:rows => self.rows, :vector_rep => @vector_rep)
        tableau.rows[row][col] = to
        return tableau

      elsif root_index == 0
        # Here we want to see if we can bubble the "to" value
        # from here, up the column to an appropriate position
        cur_row = row+1 # Start one below, as we will increment up at start
        loop do
          # TODO handle case of one row in column
          if cur_row == 0
            break
          else
            cur_row -= 1
          end

          case cur_row
          when row
            if cur_row > 0 and @rows[cur_row-1] and @rows[cur_row-1][col]
              next if not verify_up_down(@rows[cur_row-1][col], to)
            end
            if @rows[cur_row+1] and @rows[cur_row+1][col]
              next if not verify_up_down(@rows[cur_row+1][col], to)
            end
          when 0
            next if not verify_up_down(to, @rows[cur_row][col])
          when @rows.size-1
            next if not verify_up_down(@rows[cur_row][col], to)
          else
            if cur_row > 0 and @rows[cur_row-1] and @rows[cur_row-1][col]
              next if not verify_up_down(@rows[cur_row-1][col], to)
            end
            if @rows[cur_row][col]
              next if not verify_up_down(@rows[cur_row][col], to)
            end
          end

          # Success at this row!
          tableau = YoungTableau.new(:rows => self.rows, :vector_rep => @vector_rep)
          if cur_row == row
            tableau.rows[cur_row][col] = to
          else
            cur_row2 = row
            loop do
              tableau.rows[cur_row2][col] = tableau.rows[cur_row2-1][col]
              cur_row2 -= 1
              if cur_row2 == cur_row
                tableau.rows[cur_row2][col] = to
                break
              end
            end
          end

          cur_row2 = 0
          success = true
          loop do
            if not tableau.verify_row_col(cur_row2, col)
              success = false
              break
            end
            cur_row2 += 1
            break if not tableau.rows[cur_row2]
          end
          
          next if not success

          return tableau
        end
      else
        return nil
      end
    end

    def verify_row_col(row, col)
      return false if not @rows[row]
      return verify_row_col_to(row, col, @rows[row][col])
    end

    def verify_row_col_to(row, col, to)
      if @rows[row+1]
        return false if not verify_up_down(to, @rows[row+1][col])
      end
      if row > 0 and @rows[row-1]
        return false if not verify_up_down(@rows[row-1][col], to)
      end
      return false if not verify_left_right(to, @rows[row][col+1])
      return false if not verify_left_right(@rows[row][col-1], to)

      return true
    end

    def verify_bubble(col, row1, row2)
    end

    # Verify that the two vectors left and right are valid
    def verify_left_right(left, right)
      return true if left.nil? or right.nil?
      return (left <= right) ? true : false
    end

    # Verify that the two vectors up and down are valid
    def verify_up_down(up, down)
      return true if up.nil? or down.nil?

      case @vector_rep.algebra.alg
      when :alg_A, :alg_C
        return (up < down) ? true : false
      when :alg_B
        return (up <= down) ? true : false
      when :alg_D
        case (up <=> down)
        when -1
          return true
        when 1
          return false
        when 0
          return (@vector_rep.node_to_level_hash[up] == @vector_rep.algebra.rank-1) ? true : false
        else
          return true
        end
      else
        return true
      end
    end

  end
end

