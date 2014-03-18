module RubyLie
  
  class HighestWeightRep
    include Representation

    attr_reader :highest_weight
    
    def initialize(highest_weight, use_young_tableau = true)
      if not highest_weight.is_a? RubyLie::Vector
        raise TypeError, "#{highest_weight.class} is not an instance of RubyLie::Vector"
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
      generate_tree(use_young_tableau)
      
      setup()
    end
    
    # Calculated via equation dim(Lambda) = Product_{positive roots} (alpha * (Lambda + rho)) / (alpha * rho)
    def dimension
      @highest_weight.dimension
    end


    def size_w_multiplicities
      return size()
      self.inject(0) do |res, node|
        res + multiplicity(node)
      end
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
            sum += (node.weight+m*positive_root)*positive_root * multiplicity(above)
            m += 1
          else
            break
          end
        end
      end

      coeff*sum
    end

  #protected
    
    def generate_tree(use_young_tableau)
      if use_young_tableau
        top_tableau = RubyLie::YoungTableau.from_highest_weight(@highest_weight)
      else
        top_tableau = nil
      end

      @root = Node.new(:weight => @highest_weight, :representation => self, :young_tableau => top_tableau) 
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
                if node_to_check.weight == vec
                  if use_young_tableau
                    # Check if Young tableau also agree
                    if node_to_check.young_tableau == cur_node.young_tableau.next_tableau(i)
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

              # If found an overlapping node, use that instead of making a new one
              if overlapping_node
                cur_node.add_child_from_simple_root(overlapping_node, i)
                overlapping_node.add_parent_from_simple_root(cur_node, i)

                #@node_to_multiplicity[overlapping_node] = multiplicity(overlapping_node)
              else
                vec_node = Node.new(:weight => vec, :representation => self, :young_tableau => (top_tableau.nil? or not use_young_tableau) ? nil : cur_node.young_tableau.next_tableau(i)) 

                cur_node.add_child_from_simple_root(vec_node, i)
                vec_node.add_parent_from_simple_root(cur_node, i)

                @levels[cur_level+1] << vec_node

                #@node_to_multiplicity[vec_node] = multiplicity(vec_node)
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
    
  end
  
end
