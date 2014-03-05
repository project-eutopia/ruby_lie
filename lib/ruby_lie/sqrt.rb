module RubyLie
  
  class Sqrt
    include Comparable
    
    attr_reader :sqrt_part
    attr_reader :front_part
    
    def self.of(n)
      return Sqrt.new(n, 1).simplify()
    end

  # should not use publicly -- hide simplification
  protected
    # when negative, choose to take the +i, not -i part (just like we
    # implicitly choose to use the positive root)
    def initialize(sqrt_part, front_part = 1)
      if sqrt_part < 0
        sqrt_part *= -1
        front_part *= Complex(0,1)
      end
      @sqrt_part = sqrt_part
      @front_part = front_part
    end
    
  public
    def simplify
      if @sqrt_part < 0
        @sqrt_part *= -1
        @front_part *= Complex(0,1)
      end
      
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
    
    def abs
      return Sqrt.new(@sqrt_part, @front_part.abs)
    end
    
    def to_f
      return @front_part * Math.sqrt(@sqrt_part)
    end
    
    def <=>(other)
      if @sqrt_part == other
        return 0
      else
        return self.to_f <=> other.to_f
      end
    end
    
    def to_s
#      return "(#{@front_part}) * (sqrt(#{@sqrt_part}))"
      if @sqrt_part == 1
        "#{@front_part}"
      elsif @sqrt_part == 0 or @front_part == 0
        "0"
      else
        if @front_part == 1
          "sqrt(#{@sqrt_part})"
        else
          "#{@front_part} * sqrt(#{@sqrt_part})"
        end
      end
    end
    
    def to_latex
      if @sqrt_part == 1
        "#{@front_part.to_latex}"
      elsif @sqrt_part == 0 or @front_part == 0
        "0"
      else
        if @front_part == 1
          "\\sqrt{#{@sqrt_part.to_latex}}"
        else
          "#{@front_part.to_latex} \\sqrt{#{@sqrt_part.to_latex}}"
        end
      end
    end
    
    def -@
#      puts "-(#{self})"
      return Sqrt.new(@sqrt_part, -@front_part)
    end
    
    def +@
#      puts "(#{self})"
      return self.copy
    end
    
    # Just simply add as regular Numeric (no factoring and such)
    def +(a)
#      puts "(#{self}) + (#{a})"
      if self == 0
        return a
      elsif a == 0
        return self
      else
        case a
        when Numeric
          if @sqrt_part == 1
            return @front_part + a
          else
            return self.to_f + a.to_f
          end
        when Sqrt
          if @sqrt_part == a.sqrt_part
            return Sqrt.new(@sqrt_part, @front_part + a.front_part)
          else
            return self.to_f + a.to_f
          end
        else
          return self.to_f + a.to_f
        end
      end
    end
    
    def -(a)
#      puts "(#{self}) - (#{a})"
      return self + (-a)
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
    
    def /(a) #/
      if a.is_a? Integer
        return Sqrt.new(@sqrt_part, @front_part * Rational(1,a)).simplify()
      elsif a.is_a? Numeric
        return Sqrt.new(@sqrt_part, @front_part / a).simplify()
      elsif a.is_a? RubyLie::Sqrt
        return Sqrt.new(@sqrt_part*a.sqrt_part, @front_part / (a.sqrt_part*a.front_part)).simplify()
      else
        self.to_f / a.to_f
      end
    end
    
    def **(a)
      return Sqrt.new(@sqrt_part ** a, @front_part ** a).simplify()
    end

    def quo(a)
      return self / a
    end
    
    def ==(a)
      s = self.simplify()
      if a.is_a? Sqrt
        a = a.simplify()
      end
      
      case s
      when Sqrt
        case a
        when Sqrt
          return (s.sqrt_part == a.sqrt_part and s.front_part == a.front_part)
        when Numeric
          return (s.sqrt_part == 1 and s.front_part == a) 
        else
          raise TypeError, "#{a.class} (#{a}) must be Numeric, Sqrt, or Matrix"
        end
      when Numeric
        case a
        when Sqrt
          return (a.sqrt_part == 1 and a.front_part == s)
        when Numeric
          return (s == a)
        else
          raise TypeError, "#{a.class} (#{a}) must be Numeric, Sqrt, or Matrix"
        end
      end
      
      case a
      when Sqrt
        s = self.simplify()
        a = a.simplify()
        if s.sqrt_part == 0 or s.front_part == 0
          
        end
        return (s.sqrt_part == a.sqrt_part and s.front_part == a.front_part)
      when Numeric
        s = self.simplify()
        return (s.sqrt_part == 1 and s.front_part == a)
      else
        raise TypeError, "#{a.class} (#{a}) must be Numeric, Sqrt, or Matrix"
      end
    end
    
    def coerce(other)
      case other
      when Numeric
#        puts "coerce #{Sqrt.new(1,other)} into #{self.copy}"
        return Sqrt.new(1,other), self.copy
      else
        raise TypeError, "#{other.class} can't be coerced into #{self.class}"
      end
    end
  end
end  
