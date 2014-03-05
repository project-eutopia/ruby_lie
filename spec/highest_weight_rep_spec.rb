require 'spec_helper'

describe RubyLie::HighestWeightRep do
  MAX_DIMENSION = 80
  
  ALGEBRAS.each do |algebra|
    context "#{algebra.to_latex}" do

      # List of all representations that fall under our dimension limit consideration
      reps = Array.new
      (1..algebra.rank).each do |i|
        if algebra.omega(i).dimension <= MAX_DIMENSION
          reps[i] = RubyLie::HighestWeightRep.new(algebra.omega(i))
        end
      end

      it "dimension of highest weight with multiplicities agrees with Weyl dimension formula" do
        reps.each do |rep|
          next if rep.nil?
          expect(rep.size_w_multiplicities).to be == rep.dimension
        end
      end

      # Keep in mind that H_i = alpha_i^\\vee \\cdot H
      it "tr(H_i * H_j) == alpha_i^\\vee * alpha_j^\\vee for i=0,...,rank" do
        reps.each do |rep|
          next if rep.nil?

          e,f,h = rep.matrix_rep_efh
          factor = Rational(algebra.alpha_dual(0)**2,(h[0]*h[0]).trace)
            
          (0..algebra.rank).each do |i|
            (0..algebra.rank).each do |j|
              expect(factor*((h[i]*h[j]).trace)).to be == algebra.alpha_dual(i) * algebra.alpha_dual(j)
            end
          end
            
        end
      end

      it "tr(E_i * F_j) == 2/alpha_i^2 * delta_{i,j} for i=1,...,rank" do
        reps.each do |rep|
          next if rep.nil?

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
      end

      it "multiplicity agrees with cases when node is in middle of root chain" do
        pending "multiplicity... not sure how it works exactly yet"
        reps.each do |rep|
          next if rep.nil?

          rep.each do |node|
            cur_sum = 0
            parents = 0
            (1..algebra.rank).each do |i|
              if node.parents[i] and node.children[i]
                cur_sum += 1
              end
              if parents[i]
                parents += 1
              end
            end
            if parents == 0
              expect(rep.multiplicity(node)).to be == 1
            else
              if rep.multiplicity(node) != 1
                expect(parents).to be == rep.multiplicity(node)
              end
            end
            
            if cur_sum == 0
              expect(rep.multiplicity(node)).to be == 1
            else
              expect(cur_sum).to be == rep.multiplicity(node)
            end
          end

        end
      end
      
      def comm(a,b)
        return a*b-b*a
      end
      it "the matrix representations agree with Cartan matrix" do
        reps.each_with_index do |rep, index|
          next if rep.nil?

          e,f,h = rep.matrix_rep_efh()
            
          (0..algebra.rank).each do |i|
            (0..algebra.rank).each do |j|
              expect(comm(h[i], e[j])).to be == algebra.extended_cartan[j,i]*e[j]
            end
          end
          
          # Check more complicated cases
          (index..algebra.rank).each do |i|
            next if (algebra.omega(index) + algebra.omega(i)).dimension > MAX_DIMENSION
            
            rep2 = RubyLie::HighestWeightRep.new(algebra.omega(index) + algebra.omega(i))
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
