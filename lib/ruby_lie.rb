require "ruby_lie/version"
require 'matrix'

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

require "ruby_lie/root_poset"
require "ruby_lie/node"
require "ruby_lie/vector"
require "ruby_lie/algebra"
require "ruby_lie/highest_weight_rep"

