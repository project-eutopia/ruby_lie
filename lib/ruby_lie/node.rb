module RubyLie

  class Node
    include Comparable

    attr_accessor :weight
    attr_accessor :young_tableau

    attr_accessor :parents
    attr_accessor :children

    attr_accessor :representation

    def initialize(params)
      @weight   = params[:weight]
      @young_tableau = params[:young_tableau]

      @parents  = Hash.new
      @children = Hash.new

      @representation = params[:representation]
    end

    def ==(other)
      case other
      when Node
        if @weight == other.weight
          if @young_tableau and other.young_tableau
            return @young_tableau == other.young_tableau
          else
            return true
          end
        else
          return false
        end
      else
        raise TypeError, "#{other.class} must be Node"
      end
    end

    # Used to compare if this is "greater" than the other by using
    # the levels within the representation as the comparator,
    # where higher levels are "greater" than
    # TODO YoungTableau comparison
    def <=>(other)
      case other
      when Node
        return @representation.node_to_level_hash[self] <=> @representation.node_to_level_hash[other]
      else
        raise TypeError, "#{other.class} must be Node"
      end
    end

    def <=(other)
      return (self <=> other) == 1 ? false : true
    end

    def >=(other)
      return (self <=> other) == -1 ? false : true
    end

    def <(other)
      return (self <=> other) == -1 ? true : false
    end

    def >(other)
      return (self <=> other) == 1 ? true : false
    end

    def get_q(index)
      q = 0
      cur_node = self

      loop do
        break if cur_node.parents[index].nil?
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
        break if cur_node.children[index].nil?
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

    def to_latex
      return @young_tableau.to_latex if @young_tableau
      return @weight.dynkin_labels.to_s
    end

    def vec_weight_to_latex
      level = @representation.node_to_level_hash[self]

      case @weight.algebra.alg
      when :alg_A
        return "#{level+1}"

      when :alg_B
        case level <=> @weight.algebra.rank
        when -1
          return "#{level+1}"
        when 0
          return "0"
        when 1
          return "\\bar{#{2*@weight.algebra.rank+1 - level}}"
        end

      when :alg_C
        if level < @weight.algebra.rank
          return "#{level+1}"
        else
          return "\\bar{#{2*@weight.algebra.rank - level}}"
        end

      when :alg_D
        case level <=> (@weight.algebra.rank-1)
        when -1
          return "#{level+1}"
        when 0
          return "#{level+1}" if (@weight * @weight.algebra.alpha_dual(@weight.algebra.rank)) == 1
          return "\\bar{#{level+1}}" if (@weight * @weight.algebra.alpha_dual(@weight.algebra.rank)) == -1
        when 1
          return "\\bar{#{2*@weight.algebra.rank-1 - level}}"
        end

      end

      return ""
    end

  end
end
