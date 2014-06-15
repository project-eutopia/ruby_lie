require 'spec_helper'

describe RubyLie::Algebra do
  ALGEBRAS.each do |algebra|
    context "#{algebra.to_latex}" do

      it "converting between alpha vectors should retain the same value" do
        (0..algebra.rank).each do |i|
          # Convert to each type
          alpha = algebra.alpha(i)

          vectors = RubyLie::VECTOR_TYPES.map do |type|
            alpha.to_type(type)
          end

          # Verify equality works
          RubyLie::VECTOR_TYPES.each do |type1|
            RubyLie::VECTOR_TYPES.each do |type2|
              expect(alpha.to_type(type1)).to be == alpha.to_type(type2)#,
                  #"failed on i=" + i.to_s + ", type1(" + type1.to_s + ") should == type2(" + type2.to_s + ")"
            end
          end

          # Verify coefficients match going back and forth
          vectors.each do |vec|
            RubyLie::VECTOR_TYPES.each do |type|
              expect(vec.coeffs).to be == vec.to_type(type).to_type(vec.type).coeffs#,
                  #"failed on i=" + i.to_s + ", " + vec.to_s + " going to type " + type.to_s + " and back"
            end
          end
        end
      end

      it "converting between omega vectors should retain the same value" do
        # No omega_0 element
        (1..algebra.rank).each do |i|
          # Convert to each type
          omega = algebra.omega(i)

          vectors = RubyLie::VECTOR_TYPES.map do |type|
            omega.to_type(type)
          end

          # Verify equality works
          RubyLie::VECTOR_TYPES.each do |type1|
            RubyLie::VECTOR_TYPES.each do |type2|
              expect(omega.to_type(type1)).to be == omega.to_type(type2)#,
                  #"failed on i=" + i.to_s + ", type1(" + type1.to_s + ") should == type2(" + type2.to_s + ")"
            end
          end

          # Verify coefficients match going back and forth
          vectors.each do |vec|
            RubyLie::VECTOR_TYPES.each do |type|
              expect(vec.coeffs).to be == vec.to_type(type).to_type(vec.type).coeffs#,
                  #"failed on i=" + i.to_s + ", " + vec.to_s + " going to type " + type.to_s + " and back"
            end
          end
        end
      end

      it "alpha_i = Cartan_ij * omega_j" do
        (1..algebra.rank).each do |i|
          alpha_test = (1..algebra.rank).inject(0) do |res, elem|
            res + algebra.cartan[i-1,elem-1] * algebra.omega(elem)
          end

          expect(alpha_test).to be == algebra.alpha(i)
        end
      end

      it "alpha_i dot alpha_j^\\vee gives Cartan matrix" do
        (1..algebra.rank).each do |i|
          (1..algebra.rank).each do |j|
            expect(algebra.alpha(i) * algebra.alpha(j, :dual => true)).to be == algebra.cartan[i-1,j-1]#,
                #"(row,col) = (" + i.to_s + "," + j.to_s + ")"
          end
        end
      end

      it "alpha_i dot omega_j^\\vee == delta_ij" do
        (1..algebra.rank).each do |i|
          (1..algebra.rank).each do |j|
            expect(algebra.alpha(i) * algebra.omega(j, :dual => true)).to be == (i==j ? 1 : 0)#,
                #"(row,col) = (" + i.to_s + "," + j.to_s + ")"
          end
        end
      end

      it "alpha_i^\\vee dot omega_j == delta_ij" do
        (1..algebra.rank).each do |i|
          (1..algebra.rank).each do |j|
            expect(algebra.alpha(i, :dual => true) * algebra.omega(j)).to be == (i==j ? 1 : 0)#,
                #"(row,col) = (" + i.to_s + "," + j.to_s + ")"
          end
        end
      end

      it "weyl vector check" do
        weyl = algebra.weyl_vector

        # Verify sum of fundamental roots
        weyl_test = (1..algebra.rank).inject(0) do |res, elem|
          res + algebra.omega(elem)
        end
        expect(weyl_test).to be == weyl

        (1..algebra.rank).each do |i|
          expect(weyl * algebra.alpha(i, :dual => true)).to be == 1#,
              #"alpha_"+i.to_s+" dual dots to 1"
        end
        expect(weyl * algebra.alpha(0, :dual => true)).to be == (1-algebra.coxeter_number(:dual => true))#,
            #"alpha_0 dual dots to 1-h^\\vee"
      end

      it "dual weyl vector check" do
        weyl_dual = algebra.weyl_vector(:dual => true)

        # Verify sum of fundamental dual roots
        weyl_dual_test = (1..algebra.rank).inject(0) do |res, elem|
          res + algebra.omega(elem, :dual => true)
        end
        expect(weyl_dual_test).to be == weyl_dual

        (1..algebra.rank).each do |i|
          expect(weyl_dual * algebra.alpha(i)).to be == 1#,
              #"alpha_"+i.to_s+" dots to 1"
        end
        expect(weyl_dual * algebra.alpha(0)).to be == (1-algebra.coxeter_number)#,
            #"alpha_0 dual dots to 1-h"
      end

      it "check of alpha^2 values" do
        case algebra.alg
        when :alg_A
          (0..algebra.rank).each do |i|
            alpha = algebra.alpha(i)
            expect(alpha*alpha).to be == 2#,
                #"alpha_"+i.to_s
          end

        when :alg_B
          (0..algebra.rank).each do |i|
            alpha = algebra.alpha(i)
            if algebra.rank == 1 or i < algebra.rank
              expect(alpha*alpha).to be == 2#,
                  #"alpha_"+i.to_s
            else
              expect(alpha*alpha).to be == 1#,
                  #"alpha_"+i.to_s
            end
          end

        when :alg_C
          (0..algebra.rank).each do |i|
            alpha = algebra.alpha(i)
                if algebra.rank == 1 or i == 0 or i == algebra.rank
              expect(alpha*alpha).to be == 2#,
                  #"alpha_"+i.to_s
            else
              expect(alpha*alpha).to be == 1#,
                  #"alpha_"+i.to_s
            end
          end

        end
      end

      poset = algebra.root_poset
      it "weyl reflection is nilpotent" do
        orig_simple_roots = (0..algebra.rank).map {|i| algebra.alpha(i)}

        # Permute by each of the possible positive roots
        poset.each do |positive_root|
          orig_simple_roots.each do |simple|
            expect(simple).to be == simple.weyl_reflect_by(positive_root).weyl_reflect_by(positive_root)
          end
        end
      end

      it "weyl reflection by self gives minus self" do
        (0..algebra.rank).each do |i|
          simple_root = algebra.alpha(i)
          expect(simple_root).to be == -simple_root.weyl_reflect_by(simple_root)
        end
      end

      it "weyl reflection preserves inner product" do
        poset.each do |positive_root|
          (0..algebra.rank).each do |i|
            alphai = algebra.alpha(i)
            (0..algebra.rank).each do |j|
              alphaj = algebra.alpha(j)
              expect(alphai*alphaj).to be == alphai.weyl_reflect_by(positive_root) * alphaj.weyl_reflect_by(positive_root)
            end
          end
        end
      end

      poset = algebra.root_poset
      it "half of sum of positive roots == weyl vector" do
        sum_positive = 0
        poset.each do |root|
          sum_positive += root
        end
        expect(Rational(1,2) * sum_positive).to be == algebra.weyl_vector
      end


      case algebra.alg
      when :alg_A
        it "has Coxeter number r+1="+(algebra.rank+1).to_s do
          algebra.coxeter_number.should == algebra.rank+1
        end

      when :alg_B
        it "has Coxeter number 2r="+(2*algebra.rank).to_s do
          algebra.coxeter_number.should == 2*algebra.rank
        end

        if algebra.rank > 1
          it "has dual Coxeter number 2r-1="+(2*algebra.rank-1).to_s do
            algebra.coxeter_number(:dual => true).should == 2*algebra.rank - 1
          end
        end

      when :alg_C
        it "has Coxeter number 2r="+(2*algebra.rank).to_s do
          algebra.coxeter_number.should == 2*algebra.rank
        end

        if algebra.rank > 1
          it "has dual Coxeter number r+1="+(algebra.rank+1).to_s do
            algebra.coxeter_number(:dual => true).should == algebra.rank+1
          end
        end

      when :alg_D
        if algebra.rank == 2
          it "has Coxeter number 3" do
            algebra.coxeter_number.should == 3
          end
        elsif algebra.rank > 2
          it "has Coxeter number 2r-2="+(2*algebra.rank-2).to_s do
            algebra.coxeter_number.should == 2*algebra.rank-2
          end
        end

      end


    end
  end
end
