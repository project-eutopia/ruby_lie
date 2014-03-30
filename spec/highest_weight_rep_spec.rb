require 'spec_helper'

describe RubyLie::HighestWeightRep do
  MAX_DIMENSION = 180
  
  ALGEBRAS.each do |algebra|
    context "Algebra = (#{algebra.to_latex})" do

      # List of all representations that fall under our dimension limit consideration
      reps = Array.new
      (1..algebra.rank).each do |i|
        if algebra.omega(i).dimension <= MAX_DIMENSION
          reps[i] = RubyLie::HighestWeightRep.new(algebra.omega(i))
        end
      end

      reps.each_with_index do |rep, index|
        next if rep.nil?

        context "Highest weight = (#{rep.highest_weight})" do

          it "POSTULATION: eigenvalue ratios with and without sqrt(coxeter_label) weighting is the same" do
            ratio = nil
            
            eigen_w_coxeter  = rep.sum_of_coxeter_weighted_matrix_reps.to_f.eigen.eigenvalues.max {|a,b| a.real <=> b.real}
            eigen_wo_coxeter = rep.sum_of_simple_roots_matrix.to_f.eigen.eigenvalues.max {|a,b| a.real <=> b.real}
            test_ratio = eigen_w_coxeter / eigen_wo_coxeter
            
            if ratio.nil?
              ratio = test_ratio
            else
              (test_ratio).should be_within(0.0000001).of(ratio)
            end
          end

          it "dimension of highest weight with multiplicities agrees with Weyl dimension formula" do
            expect(rep.size_w_multiplicities).to be == rep.dimension
          end

          # Keep in mind that H_i = alpha_i^\\vee \\cdot H
          it "tr(H_i * H_j) == alpha_i^\\vee * alpha_j^\\vee for i=0,...,rank" do
            e,f,h = rep.matrix_rep_efh
            factor = Rational(algebra.alpha_dual(0)**2,(h[0]*h[0]).trace)
              
            (0..algebra.rank).each do |i|
              (0..algebra.rank).each do |j|
                expect(factor*((h[i]*h[j]).trace)).to be == algebra.alpha_dual(i) * algebra.alpha_dual(j)
              end
            end
                
          end

          it "tr(E_i * F_j) == 2/alpha_i^2 * delta_{i,j} for i=1,...,rank" do
            e,f,h = rep.matrix_rep_efh
            factor = Rational(algebra.alpha_dual(0)**2,(h[0]*h[0]).trace)
            
            (1..algebra.rank).each do |i|
              (1..algebra.rank).each do |j|
                if i==j
                  expect(factor*((e[i]*f[j]).trace)).to be == Rational(2,algebra.alpha(i)**2)
                else
                  expect(factor*((e[i]*f[j]).trace)).to be == 0
                end
                expect(factor*((e[i]*e[j]).trace)).to be == 0
              end
            end
                
          end

          
          def comm(a,b)
            return a*b-b*a
          end
          it "the matrix representations agree with Cartan matrix" do
            e,f,h = rep.matrix_rep_efh()
              
            (0..algebra.rank).each do |i|
              (0..algebra.rank).each do |j|
                expect(comm(h[i], e[j])).to be == algebra.extended_cartan[j,i]*e[j]
              end
            end
            
            # Check more complicated cases
            (index..algebra.rank).each do |index2|
              next
              next if (algebra.omega(index) + algebra.omega(index2)).dimension > MAX_DIMENSION
              
              rep2 = RubyLie::HighestWeightRep.new(algebra.omega(index) + algebra.omega(index2))
              e,f,h = rep2.matrix_rep_efh
                
              (1..algebra.rank).each do |i|
                (1..algebra.rank).each do |j|
                  expect(comm(h[i], e[j])).to be == algebra.extended_cartan[j,i]*e[j]
                end
              end
            end
          end
        end
      end

      
    end
  end
end
