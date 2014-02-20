module RubyLie
  
  class Sqrt
    
    attr_reader :sqrt_part
    attr_reader :front_part
    
    def self.of(n, front = 1)
      return Sqrt.new(n, front).simplify()
    end

    # TODO should not use publicly
    def initialize(sqrt_part, front_part = 1)
      @sqrt_part = sqrt_part
      @front_part = front_part
    end
    
    def simplify
      if @sqrt_part.is_a? Float
        return @front_part * Math.sqrt(@sqrt_part)
      end
      return @front_part if @sqrt_part == 1
      return 0 if @front_part == 0 or @sqrt_part == 0

      i = 2
      while i**2 <= @sqrt_part
        while @sqrt_part % (i**2) == 0
          @sqrt_part /= i**2
          @front_part *= i
          return @front_part if @sqrt_part == 1
        end
        i += 1
      end

      return self
    end
    
    def copy
      return Sqrt.new(@sqrt_part, @front_part)
    end
    
    def to_f
      return @front_part * Math.sqrt(@sqrt_part)
    end
    
    def to_s
      if @sqrt_part != 1
        if @front_part != 1
          "#{@front_part} * sqrt(#{@sqrt_part})"
        else
          "sqrt(#{@sqrt_part})"
        end
      else
        "#{@front_part}"
      end
    end 
    
    # Just simply add as regular Numeric (no factoring and such)
    def +(a)
      return self.to_f + a
    end
    
    def *(a)
      case a
      when Sqrt
        return Sqrt.new(@sqrt_part*a.sqrt_part, a.front_part*@front_part).simplify()
      when Numeric
        return Sqrt.new(@sqrt_part, a*@front_part).simplify()
      when Matrix
        return a.map do |item|
          self * item
        end
      else
        raise TypeError, "#{a.class} (#{a}) must be Numeric, Sqrt, or Matrix"
      end
    end
    
    def **(a)
      return Sqrt.new(@sqrt_part ** a, @front_part ** a).simplify()
    end
    
    def ==(a)
      case a
      when Sqrt
        return (self.sqrt_part == a.sqrt_part and self.front_part == a.front_part)
      when Numeric
        return (self.sqrt_part == 1 and self.front_part == a)
      else
        raise TypeError, "#{a.class} (#{a}) must be Numeric, Sqrt, or Matrix"
      end
    end
    
    def coerce(other)
      if @sqrt_part == 1
        return other, @front_part
      end
      
      case other
      when Numeric
        return Sqrt.new(1,other), self
      else
        raise TypeError, "#{other.class} can't be coerced into #{self.class}"
      end
    end
  end
end  
