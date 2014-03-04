module RubyLie
  
  class HighestWeightRep
    include Enumerable
    
    attr_reader :highest_weight
    attr_reader :algebra
    attr_reader :chains
    
    def initialize(highest_weight)
      if not highest_weight.is_a? RubyLie::Vector
        raise TypeError "#{highest_weight.class} is not an instance of RubyLie::Vector"
      end
      
      highest_weight = highest_weight.to_type(:omega)
      highest_weight.coeffs.each do |coeff|
        if coeff < 0
          raise HighestWeightNotDominant, "Negative coefficient #{coeff}"
        elsif not coeff.is_a? Integer
          raise HighestWeightNotDominant, "#{coeff} not an Integer" 
        end
      end
      
      @highest_weight = highest_weight.copy
      @algebra = @highest_weight.algebra
      
      # sets up @root, and @levels variables
      generate_tree
      
      # Hash from Node to index (the total index, from 0 to rank(rep)-1)
      @node_to_index_hash = Hash.new
      self.each_with_index do |node, index|
        @node_to_index_hash[node] = index
      end
      
      affinize_representation

      @chains = Hash.new
      (0..@algebra.rank).each do |i|
        @chains[i] = self.chain_hash(i)
      end
    end
    
    def size
      @node_to_index_hash.size
    end

    # Calculated via equation dim(Lambda) = Product_{positive roots} (alpha * (Lambda + rho)) / (alpha * rho)
    def dimension
      @highest_weight.dimension
    end

    def size_w_multiplicities
      self.inject(0) do |res, node|
        res + multiplicity(node)
      end
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
    def affinize_representation
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
            node_to_compare.weight == new_weight
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
      
      return Matrix.build(size, size) do |row,col|
        if hash_of_rowcol_to_val[[row,col]]
          hash_of_rowcol_to_val[[row,col]]
        else
          0
        end
      end
    end

    def get_hash_of_rowcol_to_val(i)
      size = @node_to_index_hash.size
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

    def multiplicity(node)
      return @node_to_multiplicity[node] if @node_to_multiplicity[node]
      return 1 if node.weight == @highest_weight
      
      coeff = Rational(2,(@highest_weight+@algebra.weyl_vector)**2 - (node.weight+@algebra.weyl_vector)**2)

      sum = 0
      
      @algebra.root_poset.each do |positive_root|
        m = 1
        loop do
          # Find the node corresponding to this + m*positive_root
          above = self.find do |node_to_compare|
            node_to_compare.weight == (node.weight + m*positive_root)
          end

          if above
            # Found it!
            # TODO super very inefficient!
            sum += (node.weight+m*positive_root)*positive_root * multiplicity(above)
            m += 1
          else
            break
          end
        end
      end

      coeff*sum
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
  
  #protected
    
    def generate_tree
      @root = Node.new(@highest_weight)
      @levels = Array.new
      @levels[0] = [@root]

      @node_to_multiplicity = Hash.new
      
      cur_level = 0
      loop do
        has_next_level = false
        
        # Loop through each weight at this level
        @levels[cur_level].each do |cur_node|
          
          # Check if we can decrement this weight by a given simple root
          # weight * alpha_i^\\vee = -(q + p) with q >= 0, p <= 0, and
          #   weight + p*alpha_i , ... weight + q*alpha_i
          # form a root chain
          (1..@algebra.rank).each do |i|
            # Get q for the current node and simple root
            q = cur_node.get_q(i)
            # -p = q + weight * alpha_i^\\vee
            # if -p > 0, then we can subtract this simple root
            if (q + cur_node.weight * @algebra.alpha_dual(i)) > 0
              has_next_level = true
              vec = cur_node.weight - @algebra.alpha(i)
              
              if @levels[cur_level+1].nil?
                @levels[cur_level+1] = Array.new
              end
              
              # Check if this weight is already in then next level
              overlapping_node = @levels[cur_level+1].find do |node_to_check|
                # TODO change check to make sure that this weight is arrived at
                # by the same simple roots?
                node_to_check.weight == vec
              end

              # If found an overlapping node, use that instead of making a new one
              if overlapping_node
                cur_node.add_child_from_simple_root(overlapping_node, i)
                overlapping_node.add_parent_from_simple_root(cur_node, i)

                @node_to_multiplicity[overlapping_node] = multiplicity(overlapping_node)
              else
                vec_node = Node.new(vec)

                cur_node.add_child_from_simple_root(vec_node, i)
                vec_node.add_parent_from_simple_root(cur_node, i)

                @levels[cur_level+1] << vec_node

                @node_to_multiplicity[vec_node] = multiplicity(vec_node)
              end
            end
          end
        end
        
        # If found a simple root we can subtract the weight by, continue
        if has_next_level
          cur_level += 1
        else
          break
        end
      end
    end
    
    def root
      @root
    end
    
    def levels
      @levels
    end
  end
  
end
