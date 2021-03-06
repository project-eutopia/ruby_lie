module RubyLie

# TODO: write mixin for "Representation" which deals with different types, like HighestWeightRep and RootPoset

#TODO write RootPoset to_latex, or BETTER YET, WRITE IT IN REPRESENTATION, then print it out in
#algebra document so that I can see which order of roots gives the highest weight from basic
#elements: like [E_3, [E_2, E_1]]
  class RootPoset
    include Enumerable

    attr_reader :algebra, :levels

    def initialize(algebra, simple_roots = nil)
      if algebra.is_a? RubyLie::Algebra
        @algebra = algebra

      if simple_roots and simple_roots.is_a? Array
        @simple_roots = simple_roots
      else
        @simple_roots = (1..@algebra.rank).map { |i| @algebra.alpha(i) }
      end

      # TODO F and G are fixed rank, and check rank of E type also
      elsif algebra[:alg] and algebra[:rank]
        @algebra = RubyLie::Algebra.new(algebra[:alg], algebra[:rank])
      else
        raise TypeError "#{highest_weight.class} is not an instance of RubyLie::Algebra or appropriate hash"
      end

      # sets up @root, and @levels variables
      generate_tree
    end

    def highest_root
      # TODO for now, to handle the case of D_2, add all those at the top level
      @levels[@levels.length-1].inject(0) do |res, node|
        res + node.weight
      end
    end

    def each
      (2..@levels.length-1).each do |level|
        @levels[level].each do |node|
          yield node.weight
        end
      end
    end

  protected

    def generate_tree
      @levels = Array.new
      @levels[0] = Array.new
      @levels[1] = Array.new
      @levels[2] = Array.new
      (1..@algebra.rank).each do |i|
        @levels[0][i-1] = RubyLie::Node.new(:weight => -@simple_roots[i-1], :representation => self)
        @levels[1][i-1] = RubyLie::Node.new(:weight => 0*@simple_roots[i-1], :representation => self)
        @levels[2][i-1] = RubyLie::Node.new(:weight => @simple_roots[i-1], :representation => self)
      end

      cur_level = 2
      loop do
        has_next_level = false

        # Loop through each weight at this level
        @levels[cur_level].each do |cur_node|

          # Check if we can increment this node by a given simple root
          # cur_vec * alpha_i^\\vee = -(q + p) with q >= 0, p <= 0, and
          #   weight + p*alpha_i , ... weight + q*alpha_i
          # form a root chain
          (1..@algebra.rank).each do |i|
            # Get p for the current node and simple root
            p = cur_node.get_p(i)
            # q = -p - weight * alpha_i^\\vee
            # if q > 0, then we can add this simple root
            if (-p - cur_node.weight * @simple_roots[i-1].dual) > 0
              has_next_level = true
              vec = cur_node.weight + @simple_roots[i-1]

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
                cur_node.add_parent_from_simple_root(overlapping_node, i)
                overlapping_node.add_child_from_simple_root(cur_node, i)
              else
                vec_node = RubyLie::Node.new(:weight => vec, :representation => self)

                cur_node.add_parent_from_simple_root(vec_node, i)
                vec_node.add_child_from_simple_root(cur_node, i)

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

  end

end
