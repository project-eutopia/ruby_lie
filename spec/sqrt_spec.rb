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
  
  it "<=> works" do
    expect(RubyLie::Sqrt.of(2) <=> 1*RubyLie::Sqrt.of(2)).to be == 0
    expect(RubyLie::Sqrt.of(2) <=> 2).to be == -1
    expect(RubyLie::Sqrt.of(2) <=> 1).to be == 1
    expect(RubyLie::Sqrt.of(2)*RubyLie::Sqrt.of(2) <=> 2).to be == 0
  end

  it "roots of negative numbers square to negative" do
    expect(RubyLie::Sqrt.of(-1)).to be == Complex(0,1)
    expect(RubyLie::Sqrt.of(-1)**2).to be == -1
  end
end
