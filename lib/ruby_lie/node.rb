module RubyLie
  
  class Node
    attr_accessor :weight
    attr_accessor :parents
    attr_accessor :children
    
    def initialize(weight)
      @weight   = weight
      @parents  = Hash.new
      @children = Hash.new
    end
    
    def ==(other)
      return @weight == other.weight
    end
    
    def get_q(index)
      q = 0
      cur_node = self
      
      loop do
        break if cur_node.parents.length == 0
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
    
  end
end
