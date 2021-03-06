module RubyLie

  class Algebra

    attr_reader :cartan
    attr_reader :alg
    attr_reader :rank
    attr_reader :alpha_to_ortho_matrix

    def initialize(alg, rank)
      @alg = alg
      @rank = rank

      @cartan = get_cartan
      @alpha_to_ortho_matrix = get_alpha_to_ortho_matrix
      @extended_cartan = get_extended_cartan
    end

    def vector_rep
      @vector_rep ||= RubyLie::HighestWeightRep.new(self.omega(1), false)
    end

    def root_poset
      @root_poset ||= RubyLie::RootPoset.new(self)
    end

    def to_latex
      if self.is_dual?
        return "#{self.base_type_to_char}^\\vee_#{@rank}"
      else
        return "#{self.base_type_to_char}_#{@rank}"
      end
    end

    def is_dual?
      case @alg
      when :alg_B_dual, :alg_C_dual, :alg_G_dual
        return true
      else
        return false
      end
    end

    def base_type_to_char
      case @alg
      when :alg_A
        return 'A'
      when :alg_B, :alg_B_dual
        return 'B'
      when :alg_C, :alg_C_dual
        return 'C'
      when :alg_D
        return 'D'
      when :alg_E
        return 'E'
      when :alg_F
        return 'F'
      when :alg_G, :alg_G_dual
        return 'G'
      end
    end

    def ==(a)
      return (@alg == a.alg and @rank == a.rank)
    end

    def fund_rep(i)
      return RubyLie::HighestWeightRep.new(omega(i))
    end

    def highest_root
      return -alpha(0)
    end

    # opts contains hash with :dual => true/false
    def alpha(i, opts = {})
      if i == 0
        if coxeter_label(0, opts) == 1
          return RubyLie::Vector.new(Matrix.row_vector(rank.times.map {|j| -coxeter_label(j+1, opts)}),
                                    opts[:dual] ? :alpha_dual : :alpha, self)
        else
          return RubyLie::Vector.new(Matrix.row_vector(rank.times.map {|j| Rational(-coxeter_label(j+1, opts),coxeter_label(0, opts))}),
                                    opts[:dual] ? :alpha_dual : :alpha, self)
        end
      elsif i > 0 and i <= rank
        return RubyLie::Vector.new(Matrix.row_vector(rank.times.map {|j| i == (j+1) ? 1 : 0}),
                                  opts[:dual] ? :alpha_dual : :alpha, self)
      else
        raise IndexOutOfRangeException
      end
    end
    # for convenience
    def alpha_dual(i)
      return alpha(i, :dual => true)
    end

    # opts contains hash with :dual => true/false
    def omega(i, opts = {})
      if i == 0
        # no zeroth fundamental root
        return nil
      elsif i > 0 and i <= rank
        return RubyLie::Vector.new(Matrix.row_vector(rank.times.map {|j| (i-1) == j ? 1 : 0}),
                                   opts[:dual] ? :omega_dual : :omega, self)
      else
        return nil
      end
    end
    # for convenience
    def omega_dual(i)
      return omega(i, :dual => true)
    end

    def root_poset
      return RubyLie::RootPoset.new(self)
    end

    def fundamental_rep(i)
      use_young_tableau = true

      if i >= 1 and i <= rank
        return RubyLie::HighestWeightRep.new(self.omega(i), use_young_tableau)
      else
        return nil
      end
    end

    def weyl_vector(opts = {})
      return RubyLie::Vector.new(Matrix.row_vector(rank.times.map {|j| 1}),
                                 opts[:dual] ? :omega_dual : :omega, self)
    end

    # opts contains hash with :dual => true/false
    def coxeter_number(opts = {})
      (0..@rank).inject(0) do |res, elem|
        res + coxeter_label(elem, opts)
      end
    end

    # opts contains hash with :dual => true/false
    def coxeter_label(i, opts = {})
      if i < 0 or i > @rank
        raise RubyLie::IndexOutOfRangeException.new
      end

      if @alg == :alg_A or @rank == 1
        return 1

      elsif @alg == :alg_B or @alg == :alg_B_dual
        if i <= 1
          return 1
        elsif i < @rank
          return 2
        elsif i == @rank
          if (opts[:dual] and @alg == :alg_B) or (not opts[:dual] and @alg == :alg_B_dual)
            return 1
          else
            return 2
          end
        end

      elsif @alg == :alg_C
        if opts[:dual]
          return 1
        else
          if i == 0 or i == @rank
            return 1
          else
            return 2
          end
        end

      elsif @alg == :alg_D
        if i == 0 or i == 1 or i == @rank-1 or i == @rank
          return 1
        else
          return 2
        end

      elsif @alg == :alg_E
        case @rank
        when 6
          case i
          when 0
            return 1
          when 1, 5
            return 1
          when 2, 4, 6
            return 2
          when 3
            return 3
          end

        when 7
          case i
          when 0, 1
            return 1
          when 2, 6, 7
            return 2
          when 3, 5
            return 3
          when 4
            return 4
          end

        when 8
          case i
          when 0
            return 1
          when 1, 7
            return 2
          when 2, 8
            return 3
          when 3, 6
            return 4
          when 4
            return 5
          when 5
            return 6
          end
        end

      elsif @alg == :alg_F
        case i
        when 0
          return 1
        when 1
          return 2
        when 2
          return 3
        when 3
          if opts[:dual]
            return 2
          else
            return 4
          end
        when 4
          if opts[:dual]
            return 1
          else
            return 2
          end
        end

      elsif @alg == :alg_G
        if i == 0
          return 1
        elsif i == 1
          return 2
        elsif i == 2
          if opts[:dual]
            return 1
          else
            return 3
          end
        end

      else
        # TODO other algebras
      end

      puts "Error: #{@alg} #{@rank} #{i}"
      raise RubyLie::NotSupportedException.new
    end

    def extended_cartan
      @extended_cartan
    end

  end

  # These protected methods of Algebra are used within this module, but not outside
  #protected
    class Algebra

      def vec_to_type(vec, type)

        case vec.type
        # Trivial case
        when type
          return vec

        when :alpha_dual
          return RubyLie::Vector.new(vec.coeffs.map.with_index {|elem, index| elem * self.two_over_alphai_sq(index+1)},
                                     :alpha, self).to_type(type)
        when :omega_dual
          return RubyLie::Vector.new(vec.coeffs.map.with_index {|elem, index| elem * self.two_over_alphai_sq(index+1)},
                                     :omega, self).to_type(type)

        when :alpha
          case type
          when :alpha_dual
            return RubyLie::Vector.new(vec.coeffs.map.with_index {|elem, index| Rational(elem,self.two_over_alphai_sq(index+1))},
                                       :alpha_dual, self)
          when :omega_dual
            vec_omega = vec.to_type(:omega)
            return RubyLie::Vector.new(vec_omega.coeffs.map.with_index {|elem, index| Rational(elem,self.two_over_alphai_sq(index+1))},
                                       :omega_dual, self)
          when :omega
            return RubyLie::Vector.new(vec.coeffs * @cartan,
                                       :omega, self)
          when :ortho
            return RubyLie::Vector.new(vec.coeffs * @alpha_to_ortho_matrix, :ortho, self)
          end

        when :omega
          case type
          when :omega_dual
            return RubyLie::Vector.new(vec.coeffs.map.with_index {|elem, index| Rational(elem,self.two_over_alphai_sq(index+1))},
                                       :omega_dual, self)
          when :alpha_dual
            vec_alpha = vec.to_type(:alpha)
            return RubyLie::Vector.new(vec_alpha.coeffs.map.with_index {|elem, index| Rational(elem,self.two_over_alphai_sq(index+1))},
                                       :alpha_dual, self)
          when :alpha
            return RubyLie::Vector.new(vec.coeffs * @cartan.inv, :alpha, self)
          when :ortho
            return self.vec_to_type(self.vec_to_type(vec, :alpha), :ortho)
          end

        when :ortho
          # First go to alpha
          if @alg == :alg_A or @rank == 1 or @alg == :alg_G
            # In orthonormal basis, sum of coefficients must be zero
            sum = vec.coeffs.inject(0) {|res, elem| res+elem}
            # Remove redundant last element
            coeffs = vec.coeffs.select.with_index {|elem, index| index < @rank}

            # If nonzero sum, need to subtract it from each element
            if sum != 0
              coeffs = coeffs.map do |c|
                case c.class
                when Integer
                  c - Rational(sum,vec.coeffs.length)
                else
                  c - sum / vec.coeffs.length
                end
              end
            end

            coeffs_row_vector = Matrix.row_vector(coeffs)

            # Now sum of elements is zero, so can ignore the @rank+1 element
            # and now can invert the @alpha_to_ortho_matrix
            short_alpha_to_ortho_matrix = Matrix.build(@rank,@rank) do |i,j|
              @alpha_to_ortho_matrix[i,j]
            end

            return RubyLie::Vector.new(coeffs_row_vector * short_alpha_to_ortho_matrix.inv, :alpha, self).to_type(type)
          else
            return RubyLie::Vector.new(vec.coeffs * @alpha_to_ortho_matrix.inv, :alpha, self).to_type(type)
          end
        end
      end


      def inner_prod(vec1, vec2)
        case vec1.type
        when :alpha
          result = vec1.coeffs * vec2.to_type(:omega_dual).coeffs.t
        when :alpha_dual
          result = vec1.coeffs * vec2.to_type(:omega).coeffs.t
        when :omega
          result = vec1.coeffs * vec2.to_type(:alpha_dual).coeffs.t
        when :omega_dual
          result = vec1.coeffs * vec2.to_type(:alpha).coeffs.t
        else
          result = vec1.to_type(:alpha).coeffs * vec2.to_type(:omega_dual).coeffs.t
        end

        # Convert to regular number
        if result.column_size == 1 and result.row_size == 1
          return result[0,0]
        else
          return nil
        end
      end

      def two_over_alphai_sq(i)
        if i < 0 or i > @rank
          raise IndexOutOfRangeException.new
        end

        if @alg == :alg_A or @rank == 1
          return 1

        elsif @alg == :alg_B
          if i < @rank
            return 1
          else
            return 2
          end

        elsif @alg == :alg_B_dual
          if i == @rank
            return 1
          else
            return 2
          end

        elsif @alg == :alg_C
          if i == 0 or i == @rank
            return 1
          else
            return 2
          end

        elsif @alg == :alg_D
          return 1

        elsif @alg == :alg_E
          return 1

        elsif @alg == :alg_F
          if i <= 2
            return 1
          else
            return 2
          end


        elsif @alg == :alg_G
          if i == 0 or i == 1
            return 1
          else
            return 3
          end

        else
          # TODO other cases
          return 1
        end
      end

      def get_extended_cartan
        return Matrix.build(@rank+1, @rank+1) do |row, col|
          alpha(row) * alpha_dual(col)
        end
      end

      def get_cartan

        if @alg == :alg_A or @rank == 1
          return Matrix.build(@rank, @rank) do |row, col|
            if row == col
              2
            elsif row == col+1 or row == col-1
              -1
            else
              0
            end
          end

        elsif @alg == :alg_B
          return Matrix.build(@rank, @rank) do |row, col|
            if row == col
              2
            elsif [row+1,col+1] == [@rank-1, @rank]
              -2
            elsif row == col+1 or row == col-1
              -1
            else
              0
            end
          end

        elsif @alg == :alg_C or @alg == :alg_B_dual
          return Matrix.build(@rank, @rank) do |row, col|
            if row == col
              2
            elsif [row+1,col+1] == [@rank, @rank-1]
              -2
            elsif row == col+1 or row == col-1
              -1
            else
              0
            end
          end

        elsif @alg == :alg_D
          return Matrix.build(@rank, @rank) do |row, col|
            if row == col
              2
            elsif [row+1,col+1] == [@rank, @rank-1] or [row+1,col+1] == [@rank-1, @rank]
              0
            elsif [row+1,col+1] == [@rank, @rank-2] or [row+1,col+1] == [@rank-2, @rank]
              -1
            elsif row == col+1 or row == col-1
              -1
            else
              0
            end
          end

        elsif @alg == :alg_E
          return Matrix.build(@rank, @rank) do |row, col|
            if row == col
              2
            elsif row+1 == @rank
              if col == @rank-4
                -1
              else
                0
              end
            elsif col+1 == @rank
              if row == @rank-4
                -1
              else
                0
              end
            elsif row == col+1 or row == col-1
              -1
            else
              0
            end
          end

        elsif @alg == :alg_F
          return Matrix.build(@rank, @rank) do |row, col|
            if row == col
              2
            elsif row+1 == 2 and col+1 == 3
              -2
            elsif row == col+1 or row == col-1
              -1
            else
              0
            end
          end

        elsif @alg == :alg_G
          return Matrix.build(2,2) do |row,col|
            if row == col
              2
            elsif row == 0 and col == 1
              -3
            else
              -1
            end
          end

        else
          # TODO other cases
        end

        return nil
      end

      def get_alpha_to_ortho_matrix

        if @alg == :alg_A or @rank == 1
          return Matrix.build(@rank, @rank+1) do |row, col|
            if row == col
              1
            elsif row+1 == col
              -1
            else
              0
            end
          end

        elsif @alg == :alg_B
          return Matrix.build(@rank, @rank) do |row, col|
            if row == col
              1
            elsif row+1 == col
              -1
            else
              0
            end
          end

        elsif @alg == :alg_C or @alg == :alg_B_dual
          # TODO should actually scale by 1/sqrt(2)
          return Matrix.build(@rank, @rank) do |row, col|
            if row == col and row == @rank-1
              RubyLie::Sqrt.of(2)
            elsif row == col
              RubyLie::Sqrt.of(2)*Rational(1,2)
            elsif row+1 == col
              RubyLie::Sqrt.of(2)*Rational(-1,2)
            else
              0
            end
          end

        elsif @alg == :alg_D
          return Matrix.build(@rank, @rank) do |row, col|
            if row+1 == @rank
              if col+1 == @rank or col+1 == @rank-1
                1
              else
                0
              end
            elsif row == col
              1
            elsif row+1 == col
              -1
            else
              0
            end
          end

        elsif @alg == :alg_E
          return Matrix.build(@rank, @rank) do |row, col|
            if row == @rank-2
              if col != @rank-1
                Rational(-1,2)
              else
                case @rank
                when 6
                  RubyLie::Sqrt.of(3)*Rational(1, 2)
                when 7
                  RubyLie::Sqrt.of(2)*Rational(1, 2)
                when 8
                  Rational(-1,2)
                end
              end
            elsif row == @rank-1
              if col == @rank-2
                -1
              elsif col == @rank-3
                1
              else
                0
              end
            else
              if row == col
                1
              elsif row+1 == col
                -1
              else
                0
              end
            end
          end

        elsif @alg == :alg_F
          Matrix[[0, 1,-1, 0],
                 [0, 0, 1,-1],
                 [0, 0, 0, 1],
                 [Rational(1,2), Rational(-1,2), Rational(-1,2), Rational(-1,2)]]

        elsif @alg == :alg_G
          return Matrix.build(2, 3) do |row, col|
            if row == col and row == 0
              1
            elsif row+col % 2 == 1
              -1
            elsif row == 1 and col == 1
              2
            else
              0
            end
          end

        else
          # TODO other cases
        end

      end
    end

end
