require 'bundler/setup'
Bundler.setup

require 'rspec'
require File.expand_path('../../lib/ruby_lie', __FILE__)

# This file was generated by the `rspec --init` command. Conventionally, all
# specs live under a `spec` directory, which RSpec adds to the `$LOAD_PATH`.
# Require this file using `require "spec_helper"` to ensure that it is only
# loaded once.
#
# See http://rubydoc.info/gems/rspec-core/RSpec/Core/Configuration
RSpec.configure do |config|
  config.treat_symbols_as_metadata_keys_with_true_values = true
  config.run_all_when_everything_filtered = true
  config.filter_run :focus

  # Run specs in random order to surface order dependencies. If you find an
  # order dependency and want to debug it, you can fix the order by providing
  # the seed, which is printed after each run.
  #     --seed 1234
  config.order = 'random'
end



# Setup list of algebras to check in our tests
ALGEBRAS = Array.new

alg_list = [:alg_D]#, :alg_B, :alg_C, :alg_D]#, :alg_B]#, :alg_C, :alg_D, :alg_E, :alg_F, :alg_G]
ranks = [1,2,3,4]#,5,6]

alg_list.each do |alg|
  ranks.each do |rank|
    if (alg == :alg_G and rank != 2)
      next
    elsif (alg == :alg_F and rank != 4)
      next
    elsif (alg == :alg_D and rank == 2)
      next
    elsif alg == :alg_E and (rank < 6 or rank > 8)
      next
    # Only do moderate small cases for the infinite series
    elsif alg != :alg_E and (rank > 4)
      next
    end
    ALGEBRAS << RubyLie::Algebra.new(alg, rank)
  end
end
