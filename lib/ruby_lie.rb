require "ruby_lie/version"
require 'matrix'

class Numeric
  def to_latex
    self.to_s
  end
end

class String
  def to_latex
    self.to_s
  end
end

class Rational
  def to_latex
    if self.denominator == 1
      "#{self.numerator}"
    else
      "\\frac{#{self.numerator.to_latex}}{#{self.denominator.to_latex}}"
    end
  end
end

class Matrix
  def to_latex(type = :pmatrix)
    latex = "\\begin{#{type}} \n"
    self.each_with_index do |e, row, col|
      if e.is_a? RubyLie::Sqrt
        e = e.to_latex
      end

      if col+1 == self.column_size
        if row+1 != self.row_size
          latex += "#{e} \\\\\ \n"
        else
          latex += "#{e} \n"
        end
      else
        latex += "#{e} & "
      end        
    end
    latex += "\\end{#{type}}"
    return latex
  end

  def to_ortho_string(per_line = 7)
    s = nil
    count = 0
    self.each_with_index do |e, row, col|
      if e == 1
        e = ""
      elsif e == 0
        next
      end

      count += 1
      if count % per_line == 0
        s += "\\\\ \n"
      end
       
      if s.nil?
        s = "#{e.to_latex} E_{#{row+1},#{col+1}}"
      else
        s += " + #{e.to_latex} E_{#{row+1},#{col+1}}"
      end
    end
    return s
  end

  #def to_s
  #  elem_s = Array.new
  #end

  def to_f
    self.map {|e| e.to_f}
  end
end

module RubyLie
  # Algebra symbols
  ALGEBRAS = [:alg_A, :alg_B, :alg_C, :alg_D, :alg_E, :alg_F, :alg_G,
              :alg_B_dual, :alg_C_dual, :alg_G_dual]
  
  VECTOR_TYPES = [:alpha, :alpha_dual, :omega, :omega_dual, :ortho]
  
  class IndexOutOfRangeException < Exception
  end
  class NotSupportedException < Exception
  end
  class HighestWeightNotDominant < Exception
  end
end

require "ruby_lie/sqrt"
require "ruby_lie/root_poset"
require "ruby_lie/node"
require "ruby_lie/vector"
require "ruby_lie/algebra"
require "ruby_lie/representation"
require "ruby_lie/highest_weight_rep"
require "ruby_lie/young_tableau"
