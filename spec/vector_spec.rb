require 'spec_helper'

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
