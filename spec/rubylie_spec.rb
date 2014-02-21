require 'spec_helper'

describe RubyLie::Sqrt do
  it "non null" do
    expect(RubyLie::Sqrt.of(0)).not_to be_nil
    expect(RubyLie::Sqrt.of(1)).not_to be_nil
    expect(RubyLie::Sqrt.of(2)).not_to be_nil
    expect(RubyLie::Sqrt.of(3)).not_to be_nil
    expect(RubyLie::Sqrt.of(4)).not_to be_nil
  end
  
  it "square roots simplify" do
    expect(RubyLie::Sqrt.of(0)).to be == 0
    expect(RubyLie::Sqrt.of(1)).to be == 1
    expect(RubyLie::Sqrt.of(4)).to be == 2
    expect(RubyLie::Sqrt.of(9)).to be == 3
    expect(RubyLie::Sqrt.of(16)).to be == 4
  end
  
  it "can multiply by Numeric" do
    expect(1*RubyLie::Sqrt.of(1)).to be == 1
    expect(RubyLie::Sqrt.of(1)*1).to be == 1
    expect(1*RubyLie::Sqrt.of(4)).to be == 2
    expect(RubyLie::Sqrt.of(4)*1).to be == 2
    expect(2*RubyLie::Sqrt.of(1)).to be == 2
    expect(RubyLie::Sqrt.of(1)*2).to be == 2
    expect(2*RubyLie::Sqrt.of(2)).to be == RubyLie::Sqrt.of(8)
    expect(RubyLie::Sqrt.of(2)*2).to be == RubyLie::Sqrt.of(8)
    expect(0*RubyLie::Sqrt.of(3)).to be == 0
    expect(RubyLie::Sqrt.of(3)*0).to be == 0
    expect(1*RubyLie::Sqrt.of(1)).to be == 1
  end
  
  it "perfect squares simplify" do
    expect(RubyLie::Sqrt.of(0)*RubyLie::Sqrt.of(0)).to be == 0
    expect(RubyLie::Sqrt.of(1)*RubyLie::Sqrt.of(1)).to be == 1
    expect(RubyLie::Sqrt.of(2)*RubyLie::Sqrt.of(2)).to be == 2
    expect(RubyLie::Sqrt.of(3)*RubyLie::Sqrt.of(3)).to be == 3
    expect(RubyLie::Sqrt.of(4)*RubyLie::Sqrt.of(4)).to be == 4
  end
end

describe RubyLie::Vector do
  a3 = RubyLie::Algebra.new(:alg_A, 3)
  
  it "can initialize a Vector from an algebra instance or from hash" do
    expect(RubyLie::Vector.new([1], :alpha, a3)).not_to be_nil
    expect(RubyLie::Vector.new(Matrix[[1]], :alpha, a3)).not_to be_nil
    expect(RubyLie::Vector.new([1], :alpha, {alg: :alg_A, rank: 1})).not_to be_nil
    expect(RubyLie::Vector.new(Matrix[[1]], :alpha, {alg: :alg_A, rank: 1})).not_to be_nil
  end

  it "adding or subtracting 0 to alpha_i should make no difference" do
    (0..a3.rank).each do |i|
      expect(a3.alpha(i)).to be == (0 + a3.alpha(i))
      expect(a3.alpha(i)).to be == (a3.alpha(i) + 0)
      expect(a3.alpha(i)).to be == (a3.alpha(i) - 0)
    end
  end   

  it "minus of vector is -1*vector" do
    (0..a3.rank).each do |i|
      expect(-a3.alpha(i)).to be == (-1 * a3.alpha(i))
    end
  end   

  it "multiply by 0 gives 0" do
    (0..a3.rank).each do |i|
      expect(0).to be == (0 * a3.alpha(i))
      expect(0).to be == (a3.alpha(i) * 0)
    end
  end

  it "multiplication by 1 gives same result" do
    (0..a3.rank).each do |i|
      expect(1*a3.alpha(i)).to be == (a3.alpha(i))
      expect(2*a3.alpha(i)*Rational(1,2)).to be == (a3.alpha(i))
    end
  end

  it "check that 2*x - x == x" do
    (0..a3.rank).each do |i|
      expect(2*a3.alpha(i) - a3.alpha(i)).to be == (a3.alpha(i))
    end
  end
  
end

describe RubyLie::Algebra do
  
  algebra_list = [:alg_A, :alg_B, :alg_C, :alg_D]
  
  algebra_list.each do |alg|
    
    context alg.to_s do
      ranks = [1,2,3,4]
      
      ranks.each do |rank|
        algebra = RubyLie::Algebra.new(alg, rank)
        
        context ("r=" + rank.to_s) do
          
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
            weyl_test = (1..rank).inject(0) do |res, elem|
              res + algebra.omega(elem)
            end
            expect(weyl_test).to be == weyl
            
            (1..rank).each do |i|
              expect(weyl * algebra.alpha(i, :dual => true)).to be == 1#,
                  #"alpha_"+i.to_s+" dual dots to 1"
            end
            expect(weyl * algebra.alpha(0, :dual => true)).to be == (1-algebra.coxeter_number(:dual => true))#,
                #"alpha_0 dual dots to 1-h^\\vee"
          end

          it "dual weyl vector check" do
            weyl_dual = algebra.weyl_vector(:dual => true)
            
            # Verify sum of fundamental dual roots
            weyl_dual_test = (1..rank).inject(0) do |res, elem|
              res + algebra.omega(elem, :dual => true)
            end
            expect(weyl_dual_test).to be == weyl_dual

            (1..rank).each do |i|
              expect(weyl_dual * algebra.alpha(i)).to be == 1#,
                  #"alpha_"+i.to_s+" dots to 1"
            end
            expect(weyl_dual * algebra.alpha(0)).to be == (1-algebra.coxeter_number)#,
                #"alpha_0 dual dots to 1-h"
          end

          it "check of alpha^2 values" do
            case alg
            when :alg_A
              (0..rank).each do |i|
                alpha = algebra.alpha(i)
                expect(alpha*alpha).to be == 2#,
                    #"alpha_"+i.to_s
              end

            when :alg_B
              (0..rank).each do |i|
                alpha = algebra.alpha(i)
                if rank == 1 or i < rank
                  expect(alpha*alpha).to be == 2#,
                      #"alpha_"+i.to_s
                else
                  expect(alpha*alpha).to be == 1#,
                      #"alpha_"+i.to_s
                end
              end

            when :alg_C
              (0..rank).each do |i|
                alpha = algebra.alpha(i)
                    if rank == 1 or i == 0 or i == rank
                  expect(alpha*alpha).to be == 2#,
                      #"alpha_"+i.to_s
                else
                  expect(alpha*alpha).to be == 1#,
                      #"alpha_"+i.to_s
                end
              end

            end
          end
          
          poset = RubyLie::RootPoset.new(algebra)
          it "half of sum of positive roots == weyl vector" do
            sum_positive = 0
            poset.each do |root|
              sum_positive += root
            end
            expect(Rational(1,2) * sum_positive).to be == algebra.weyl_vector
          end
          
          it "highest root from poset agrees with -alpha_0" do
            expect(poset.highest_root).to be == -algebra.alpha(0)
          end
          
          def comm(a,b)
            return a*b-b*a
          end
          it "the matrix representations agree with Cartan matrix" do
            (1..algebra.rank).each do |omega_i|
              rep = RubyLie::HighestWeightRep.new(algebra.omega(omega_i))
              e,f,h = rep.matrix_rep_efh
              
              (1..algebra.rank).each do |i|
                (1..algebra.rank).each do |j|
                  expect(comm(h[i], e[j])).to be == algebra.cartan[j-1,i-1]*e[j]
                end
              end
            end
          end
          
          case alg
          when :alg_A
            it "has Coxeter number r+1="+(rank+1).to_s do
              algebra.coxeter_number.should == rank+1
            end

          when :alg_B
            it "has Coxeter number 2r="+(2*rank).to_s do
              algebra.coxeter_number.should == 2*rank
            end
            
            if rank > 1
              it "has dual Coxeter number 2r-1="+(2*rank-1).to_s do
                algebra.coxeter_number(:dual => true).should == 2*rank - 1
              end
            end

          when :alg_C
            it "has Coxeter number 2r="+(2*rank).to_s do
              algebra.coxeter_number.should == 2*rank
            end
            
            if rank > 1
              it "has dual Coxeter number r+1="+(rank+1).to_s do
                algebra.coxeter_number(:dual => true).should == rank+1
              end
            end

          when :alg_D
            if rank == 2
              it "has Coxeter number 3" do
                algebra.coxeter_number.should == 3
              end
            elsif rank > 2
              it "has Coxeter number 2r-2="+(2*rank-2).to_s do
                algebra.coxeter_number.should == 2*rank-2
              end
            end

          end
        end
      end #end ranks
      
    end #end algebra

  end
end