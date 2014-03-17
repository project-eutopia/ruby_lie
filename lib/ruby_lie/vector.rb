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
      if algebra.is_a? RubyLie::Algebra
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

    def dominant?
      dynkin_labels.each do |label|
        if label < 0
          return false
        end
      end

      return true
    end

    # Calculated via equation dim(Lambda) = Product_{positive roots} (alpha * (Lambda + rho)) / (alpha * rho)
    def dimension
      if not dominant?
        return 0
      end
      
      prod = 1
      
      poset = RubyLie::RootPoset.new(@algebra)
      rho = @algebra.weyl_vector
      
      poset.each do |root|
        prod *= Rational(root * (self + rho), root * rho)
      end

      return prod.denominator == 1 ? prod.numerator : prod
    end
    
    def type_to_latex(type)
      case type
      when :alpha
        "\\alpha"
      when :alpha_dual
        "\\alpha^\\vee"
      when :omega
        "\\omega"
      when :omega_dual
        "\\omega^\\vee"
      when :ortho
        "e"
      end
    end
    
    def to_latex
      latex = ""
      first = true
      self.coeffs.each_with_index do |label, row, col|
        if label != 0
          if first
            first = false
          else
            if label > 0
              latex += " + "
            else
              latex += " - "
            end
          end

          if label.is_a? Rational and label.denominator == 1
            label = label.numerator
          end
          
          if label.abs == 1
            label = ""
          else
            label = label.abs
          end
          latex += "#{label.to_latex} #{type_to_latex(@type)}_{#{col+1}}"
        end
      end
      
      return latex
    end

    def dual
      return (2 / (self**2)) * self
    end

    def weyl_reflect_by(alpha)
      self - (2 * (self*alpha) / alpha**2) * alpha
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

    def **(n)
      case n
      when Integer
        if n > 0
          return n.times.inject(1) {|res, elem| res * self}
        end
      end

      raise TypeError, "#{n} should be positive integer to exponentiate on #{self.class}"
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

    def zero?
      @coeffs.each do |i|
        if i != 0
          return false
        end
      end

      true
    end
    
    def ==(v)
      return false if @algebra != v.algebra
      
      v = v.to_type(@type) if @type != v.type
      
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
