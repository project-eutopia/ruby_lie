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
        yield rows[row][i], row if rows[row][i].nil?
        row -= 1
        break if row == 0
      end
    end


    def ==(tableau)
      case tableau
      when YoungTableau
        return self.rows == tableau.rows
      when Array
        return self.rows == tableau
      else
        throw TypeError, "#{tableau.class} must be of type YoungTableau or Array"
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

    def try_col_with_root_index(col, root_index)

      if root_index >= 0
        self.each_col_at(col) do |elem, row|

          @vector_rep.each_chain_link(root_index) do |from, to|
            if from == elem and try_next_tableau(row, col, to, root_index)
              tableau = YoungTableau.new(:rows => self.rows, :vector_rep => @vector_rep)
              tableau.rows[row][col] = to
              return tableau
            end
          end
              
        end

      elsif root_index == 0
        self.each_col_at_backwards(col) do |elem, row|
          return nil
        end
      else
        return nil
      end

      return nil

    end


    def try_next_tableau(row, col, to, root_index)
      # Try this possible change
      if @rows[row][col+1] and @rows[row][col+1] < to
        return false
      elsif @rows[row+1] and @rows[row+1][col]

        # A and C type must be strictly increasing down a column
        if (@vector_rep.algebra.alg == :alg_A or @vector_rep.algebra.alg == :alg_C) and @rows[row+1][col] <= to
          return false
        elsif @vector_rep.algebra.alg == :alg_B and @rows[row+1][col] < to
          return false
        elsif @vector_rep.algebra.alg == :alg_D
          if @rows[row+1][col] < to
            return false
          # Allow the pair to be n and \bar{n}, but not other combinations
          elsif (@rows[row+1][col] <=> to) == 0 and @vector_rep.node_to_level_hash[to] != @vector_rep.algebra.rank-1
            return false
          end
        end

      end

      return true
    end
  end
end

