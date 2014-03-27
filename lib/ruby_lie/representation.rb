module RubyLie
  
  module Representation
    include Enumerable
    
    attr_reader :algebra
    attr_reader :chains
    attr_reader :levels
    attr_reader :node_to_index_hash
    attr_reader :node_to_level_hash

    # User of this class must initialize the levels and algebra
    def setup(use_young_tableau)
      # Hash from Node to index (the total index, from 0 to rank(rep)-1)
      @node_to_index_hash = Hash.new
      self.each_with_index do |node, index|
        @node_to_index_hash[node] = index
      end

      @node_to_level_hash = Hash.new
      @levels.each_with_index do |level, index|
        level.each do |node|
          @node_to_level_hash[node] = index
        end
      end

      
      affinize_representation(use_young_tableau)

      @chains = Hash.new
      (0..@algebra.rank).each do |i|
        @chains[i] = self.chain_hash(i)
      end
    end

    def size
      @node_to_index_hash.size
    end


    def node_at_index(index)
      @node_to_index_hash.each do |cur_node, cur_index|
        if index == cur_index
          return cur_node
        end
      end
      return nil
    end
    
    # Add links in the chain for \alpha_0 root
    def affinize_representation(use_young_tableau = true)
      # TODO
      # for now treat the disjoint D_2 case specially... somehow
      if @algebra.alg == :alg_D and @algebra.rank == 2
        @levels[0][0].add_child_from_simple_root(@levels[0][0], 0)
        @levels[0][0].add_parent_from_simple_root(@levels[0][0], 0)
        @levels[0][1].add_child_from_simple_root(@levels[0][1], 0)
        @levels[0][1].add_parent_from_simple_root(@levels[0][1], 0)
      end
        
      
      # Loop through all weights and find those where we can subtract alpha_0
      #
      # weight * alpha_i^\\vee = -(q + p) with q >= 0, p <= 0, and
      #   weight + p*alpha_i , ... weight + q*alpha_i
      # -p = q + weight * alpha_0^\\vee
      #
      # Going down the chain in order means that we can assume q=0 for each we find
      alpha_0 = @algebra.alpha(0)
      alpha_0_dual = @algebra.alpha(0, :dual => true)
      
      self.each do |node|
        abs_p = node.get_q(0) + node.weight * alpha_0_dual
#        puts "abs_p = #{abs_p}"

        cur_node = node
        while abs_p > 0
          new_weight = cur_node.weight - alpha_0
          
          # Find the next node
          next_node = self.find do |node_to_compare|
            #node_to_compare.weight == new_weight
            if node_to_compare.weight == new_weight
              if use_young_tableau
                # Check if Young tableau also agree
                if node_to_compare.young_tableau == cur_node.young_tableau.next_tableau(0)
                  true
                else
                  false
                end
              else
                true
              end
            else
              false
            end
          end
          
          if not next_node
            puts "#{cur_node} #{next_node}"
            break
          end

          cur_node.add_child_from_simple_root(next_node, 0)
          next_node.add_parent_from_simple_root(cur_node, 0)
          cur_node = next_node
          abs_p -= 1
        end
      end
    end
    
    def each
      @levels.each do |level|
        level.each do |cur_node|
          yield cur_node
        end
      end
    end

    def node_to_index(node)
      return @node_to_index_hash[node]
    end
    
    def each_with_index
      index = 0
      @levels.each do |level|
        level.each do |cur_node|
          yield cur_node, index
          index += 1
        end
      end
    end
    
    # Returns pairs of nodes [from, to], which denote each link on the chains
    def each_chain_link(root_index)
      (0..self.num_chains(root_index)).each do |chain_num|
        self.each_chain(root_index, chain_num) do |node|
          if node.children[root_index]
            yield node, node.children[root_index]
          end
        end
      end
    end

    def each_chain(root_index, chain_num)
      if @chains[root_index]
        @chains[root_index][:node_to_chain].each do |cur_node, cur_chain_num|
          if chain_num == cur_chain_num
            yield cur_node
          end
        end
      end
    end
    
    def each_chain_with_index(root_index, chain_num)
      index = 0
      if @chains[root_index]
        @chains[root_index][:node_to_chain].each do |cur_node, cur_chain_num|
          if chain_num == cur_chain_num
            yield cur_node, index
            index += 1
          end
        end
      end
    end
    
    # i is of alpha_i, and this gives an iterator down the whole chain with the given node (defaults to root)
    def each_in_root_chain(i, node = @root)
      if node.nil?
        raise ArgumentError, "node must be non-nil"
      elsif not node.is_a? Node
        raise ArgumentError, "node must be instance of Node"
      end
      
      # Find the top node
      loop do
        if node.parents[i]
          node = node.parents[i]
        else
          break
        end
      end
      
      # Now pointing at top, loop down, if chain exists
      if node.children[i]
        loop do
          yield node
          
          break if not node.children[i]
          node = node.children[i]
        end
      end
    end
    
    def E(i,j,n)
      Matrix.build(n,n) do |row,col|
        if row == i and j == col
          1
        else
          0
        end
      end
    end
    
    def num_chains(i)
      @chains[i][:num_chains]
    end
    
    def chain_length(root_index, chain_num)
      l = 0
      self.each_chain(root_index, chain_num) do |node|
        l += 1
      end
      return l
    end
    
    # Returns object like the following
    #
    # chain_hash[:num_chains] = number of chains
    # chain_hash[:node_to_chain] = {node1 => node1_chain_num, node2 => ...}
    def chain_hash(i)
      ret_hash = {:num_chains => 0, :node_to_chain => Hash.new}
      
      found = false
      
      self.each do |node|
        # If already added, ignore
        if not ret_hash[:node_to_chain][node]
          found = false
          self.each_in_root_chain(i, node) do |node_in_chain|
            ret_hash[:node_to_chain][node_in_chain] = ret_hash[:num_chains]
            found = true
          end
          
          # If this node was indeed part of a chain, add it
          if found
            ret_hash[:num_chains] += 1
          end
        end
      end
      
      return ret_hash
    end

    def matrix_rep(i)
      
      hash_of_rowcol_to_val = self.get_hash_of_rowcol_to_val(i)
      this_size = self.size
      
      return Matrix.build(this_size, this_size) do |row,col|
        if hash_of_rowcol_to_val[[row,col]]
          hash_of_rowcol_to_val[[row,col]]
        else
          0
        end
      end
    end

    def get_hash_of_rowcol_to_val(i)
      size = self.size
      h = Hash.new

      (0..self.num_chains(i)).each do |chain_num|
        self.each_chain_with_index(i, chain_num) do |node, node_num|
          # Skip the first node of the chain, to get (chain.length-1) nodes
          if node.parents[i]
            # This subchain is a su(2) subalgebra with spin j = (chain_length-1)/2
            # So we use the normalization:
            #   J^- |j,m> = sqrt(j(j+1) - m(m-1)) |j,m-1>
            # and note that J^- = transpose(J^+) to get the following coefficient
            h[[@node_to_index_hash[node.parents[i]],@node_to_index_hash[node]]] = Sqrt.of(node_num * (chain_length(i,chain_num) - node_num))
          end
        end
      end

      return h
    end
    
    def matrix_rep_efh
      e = Hash.new
      f = Hash.new
      h = Hash.new
      
      (0..@algebra.rank).each do |i|
        e[i] = matrix_rep(i)
        f[i] = e[i].t
        h[i] = e[i]*f[i] - f[i]*e[i]
      end
      
      return [e, f, h]
    end
    
    def sum_of_simple_roots_matrix
      res = Matrix.zero(self.size, self.size)
      (0..@algebra.rank).each do |i|
        res += matrix_rep(i)
      end
      return res
    end
    
    def sum_of_coxeter_weighted_matrix_reps
      res = Matrix.zero(self.size, self.size)
      (0..@algebra.rank).each do |i|
        res += RubyLie::Sqrt.of(@algebra.coxeter_label(i, :dual => true)) * matrix_rep(i)
      end
      return res
    end

    
    
    def to_s
      return "\#<RubyLie::HighestWeightRep \n" +
          @levels.each_with_index.inject("") {|res, (value, index)|
            res + "LEVEL #{index} :: \n" + value.inject("") {|res2, value|
              res2 + "    " + value.to_s + "\n"
            } + "\n"
          } + ">"
    end
    
    def to_latex
      latex  = "\\begin{tikzpicture}[->,>=stealth',shorten >=1pt,auto,node distance=2.2cm,\n"
      latex += "                    thick,main node/.style={draw,font=\\sffamily\\large\\bfseries}]\n"
      latex += "\n"

      latex += @levels.each_with_index.inject("") do |res, (value, index)|
        res + value.each_with_index.inject("") do |res2, (node, index2)|
          # First has no position
          if (index == 0) and (index2 == 0)
            res2 + "  \\node[main node] (h_#{index}_#{index2}) {$#{node.weight.dynkin_labels.to_s}$};\n"
            
          # First of a level just below the first of the previous level
          elsif index2 == 0
            res2 + "  \\node[main node] (h_#{index}_#{index2}) [below of=h_#{index-1}_#{index2}] {$#{node.weight.dynkin_labels.to_s}$};\n"
            
          # Right of the previous index
          else
            res2 + "  \\node[main node] (h_#{index}_#{index2}) [right of=h_#{index}_#{index2-1}] {$#{node.weight.dynkin_labels.to_s}$};\n"
          end
        end
      end
      
      latex += "\n"
      latex += "  \\path[every node/.style={font=\\sffamily\\small}]\n"

      latex += @levels.each_with_index.inject("") do |res, (level, index)|
        res + level.each_with_index.inject("") do |res2, (node, index2)|
          res2 + "    (h_#{index}_#{index2})\n" + node.children.inject("") do |res3, (key, node2)|
            # TODO for now, just skip the arrows that loop back through the alpha_0 root
            if key != 0
              s = res3 + "        edge node {$\\alpha_#{key}$} (h_#{index+1}_"
              @levels[index+1].each_with_index do |node_to_check, index3|
                if node_to_check.weight == node2.weight
                  s += "#{index3})\n"
                end
              end
              s
            else
              res3
            end
          end
        end
      end
      
      latex += "  ;\n"
      latex += "\\end{tikzpicture}\n"
    end
  end
  
end
