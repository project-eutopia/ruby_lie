require 'matrix'

module RubyLie
  
  class Vector
    attr_reader :coeffs, :type, :algebra
    
    # algebra can be an instance of RubyLie::Algebra or a hash e.g. {:alg => :alg_A, :rank => 3}
    def initialize(coeffs, type, algebra)
      case coeffs
      when Matrix
        @coeffs = coeffs
      when Array
        @coeffs = Matrix.row_vector(coeffs)
      else
        return nil
      end
      
      @type = type
      if algebra.class == RubyLie::Algebra
        @algebra = algebra
      
      # TODO F and G are fixed rank, and check rank of E type also
      elsif algebra[:alg] and algebra[:rank]
        @algebra = RubyLie::Algebra.new(algebra[:alg], algebra[:rank])
      else
        return nil
      end
    end
    
    def to_type(type)
      return @algebra.vec_to_type(self, type)
    end
    
    def dynkin_labels
      vec_in_omega = self.to_type(:omega)
      vec_in_omega.coeffs.map {|i| i}.to_a[0]
    end
    
    def +(v)
      case v
      when Numeric
        if v != 0
          raise TypeError, "cannot add non-zero number to vector"
        else
          return self
        end
      end
      
      if @type != v.type
        v = v.to_type(@type)
      end
      
      # TODO implement each, map, etc to make it easier to do stuff like this
      return RubyLie::Vector.new(@coeffs.map.with_index {|elem, index| elem + v.coeffs[0,index]}, @type, @algebra)
    end
    
    def -(v)
      return self + (-1)*v
    end
    
    def -@
      return RubyLie::Vector.new(@coeffs.map {|elem| -elem}, @type, @algebra)
    end
    
    def *(v)
      case v
      when Numeric
        if v == 0
          return 0
        else
          return RubyLie::Vector.new(@coeffs.map {|elem| v * elem}, @type, @algebra)
        end
      when RubyLie::Vector
        return @algebra.inner_prod(self, v)
      else
        raise TypeError, "#{self.class} can't be multiplied by class #{v.class}"
      end
    end

    # Borrowed from Matrix source
    #
    # The coerce method provides support for Ruby type coercion.
    # This coercion mechanism is used by Ruby to handle mixed-type
    # numeric operations: it is intended to find a compatible common
    # type between the two operands of the operator.
    # See also Numeric#coerce.
    #
    def coerce(other)
      case other
      when Numeric
        return self, other
      else
        raise TypeError, "#{self.class} can't be coerced into #{other.class}"
      end
    end
    
    def ==(v)
      if @type != v.type
        v = v.to_type(@type)
      end
      
      return @coeffs == v.coeffs
    end
    
    def to_s
      return "\#<RubyLie::Vector " +
             coeffs.to_a[0].to_s + ", :type => " + @type.to_s + 
                           ", :alg => "  + @algebra.alg.to_s + 
                           ", :rank => " + @algebra.rank.to_s +
                           ">"
    end
    
    def copy
      return RubyLie::Vector.new(@coeffs.map {|elem| elem}, @type, @algebra)
    end
  end
  
end
