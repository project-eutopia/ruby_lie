require 'matrix'

module RubyLie
  
  class HighestWeightRep
    
    attr_reader :highest_weight
    attr_reader :algebra
    
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

      @chains = Hash.new
      (1..@algebra.rank).each do |i|
        @chains[i] = self.chain_hash(i)
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

    # TODO implement highest root, and sqrt factors
    def matrix_rep(i)
      node_to_index = Hash.new

      # TODO combine into bottom loop?
      self.each_with_index do |node, index|
        node_to_index[node] = index
      end
      
      e1 = Matrix.zero(node_to_index.size)

      self.each_with_index do |node, index|
        node.parents.each do |key_alpha, above_node|
          if key_alpha == i
            e1 += E(node_to_index[above_node], index, node_to_index.size)
          end
        end
      end
      
      return e1
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
      latex =  "\\documentclass{article}\n"
      latex += "\\usepackage{tikz}\n"
      latex += "\\usetikzlibrary{arrows}\n"
      latex += "\\begin{document}\n"
      latex += "\\begin{tikzpicture}[->,>=stealth',shorten >=1pt,auto,node distance=3cm,\n"
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
            s = res3 + "        edge node {$\\alpha_#{key}$} (h_#{index+1}_"
            @levels[index+1].each_with_index do |node_to_check, index3|
              if node_to_check.weight == node2.weight
                s += "#{index3})\n"
              end
            end
            s
          end
        end
      end
      
      latex += "  ;\n"
      latex += "\\end{tikzpicture}\n"
      latex += "\\end{document}\n"
    end
  
    class Node
      attr_accessor :weight
      attr_accessor :parents
      attr_accessor :children
      
      def initialize(weight)
        @weight   = weight
        @parents  = Hash.new
        @children = Hash.new
      end
      
      def get_q(index)
        q = 0
        cur_node = self
        
        loop do
          break if cur_node.parents.length == 0
          break if cur_node.parents[index] == nil
          cur_node = cur_node.parents[index]
          q += 1
        end
        
        return q
      end
      
      def get_p(index)
        p = 0
        cur_node = self
        
        loop do
          break if cur_node.children.length == 0
          break if cur_node.children[index] == nil
          cur_node = cur_node.children[index]
          p -= 1
        end
        
        return p
      end
      
      def to_s
        weight.to_s
      end
      
      def add_parent_from_simple_root(parent, root_i)
        @parents[root_i] = parent
      end
      
      def add_child_from_simple_root(child, root_i)
        @children[root_i] = child
      end
      
    end
    
  #protected
    
    def generate_tree
      @root = Node.new(@highest_weight)
      @levels = Array.new
      @levels[0] = [@root]
      
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
              else
                vec_node = Node.new(vec)

                cur_node.add_child_from_simple_root(vec_node, i)
                vec_node.add_parent_from_simple_root(cur_node, i)

                @levels[cur_level+1] << vec_node
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