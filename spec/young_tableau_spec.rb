
require 'spec_helper'

describe RubyLie::YoungTableau do
  it "initializes from partition" do
    RubyLie::YoungTableau.new(:partition => [1]).should_not be_nil
    RubyLie::YoungTableau.new(:partition => [5]).should_not be_nil
    RubyLie::YoungTableau.new(:partition => [2,1]).should_not be_nil
    RubyLie::YoungTableau.new(:partition => [3,2,1,1]).should_not be_nil
  end
  
  it "A-type simple algebra, should increment properly" do
    a4 = RubyLie::Algebra.new(:alg_A, 4)
    vec_rep = a4.vector_rep

    tableau = RubyLie::YoungTableau.from_highest_weight(a4.omega(2))
    tableau.should_not be_nil
    tableau.partition.should be == [1,1]
    tableau.rows[0][0].should be == vec_rep.levels[0][0]
    tableau.rows[1][0].should be == vec_rep.levels[1][0]
    
    tableau.next_tableau(1).should be_nil

    next_tableau = tableau.next_tableau(2)
    next_tableau.should_not be_nil
    next_tableau.rows[0][0].should be == vec_rep.levels[0][0]
    next_tableau.rows[1][0].should be == vec_rep.levels[2][0]


    tableau.next_tableau(3).should be_nil
  end

  it "check equality of simple case" do
    a4 = RubyLie::Algebra.new(:alg_A, 4)
    vec_rep = a4.vector_rep

    tableau = RubyLie::YoungTableau.from_highest_weight(a4.omega(2))
    tableau2 = RubyLie::YoungTableau.from_highest_weight(a4.omega(2))
    tableau.should be == tableau2

    tableau.next_tableau(2).should be == tableau2.next_tableau(2)
  end

  it "omega_1 rep same as vector rep" do
    a2 = RubyLie::Algebra.new(:alg_A, 2)
    vec_rep = a2.vector_rep
    rep1 = RubyLie::HighestWeightRep.new(a2.omega(1), true)

    top_vec = vec_rep.levels[0][0]
    top_1   = rep1.levels[0][0]

    
    (top_vec).should be == top_1
    (top_vec.children[1]).should be == top_1.children[1]

    bot_vec = top_vec.children[1].children[2]
    bot_1 = top_1.children[1].children[2]

    bot_vec.should be == bot_1

    bot_vec.children[0].should_not be_nil
    bot_1.children[0].should_not be_nil
  end

  it "check omega_2 of a3" do
    a3 = RubyLie::Algebra.new(:alg_A, 3)
    vec_rep = a3.vector_rep
    top_tableau = RubyLie::YoungTableau.from_highest_weight(a3.omega(2))

    top_tableau.next_tableau(1).should be_nil
    top_tableau.next_tableau(2).should_not be_nil
    top_tableau.next_tableau(3).should be_nil

    next_tableau = top_tableau.next_tableau(2)

    next_tableau.next_tableau(1).should_not be_nil
  end
  
  it "each_col_at(_backwards) works" do
    a3 = RubyLie::Algebra.new(:alg_A, 3)
    tableau = RubyLie::YoungTableau.from_highest_weight(a3.omega(3))

    i = 0
    tableau.each_col_at(0) do |elem, row|
      elem.should be == tableau.rows[row][0]
      row.should be == i
      i += 1
    end

    i = 2
    tableau.each_col_at_backwards(0) do |elem, row|
      elem.should be == tableau.rows[row][0]
      row.should be == i
      i -= 1
    end
  end

  it "check affine simple root" do
    a3 = RubyLie::Algebra.new(:alg_A, 3)
    tableau = RubyLie::YoungTableau.from_highest_weight(a3.omega(2))

    second_last = tableau.next_tableau(2).next_tableau(1).next_tableau(3)
    second_last.should_not be_nil
    second_last.next_tableau(0).should be == tableau
    second_last.next_tableau(2).next_tableau(0).should be == tableau.next_tableau(2)
  end

  it "spinor reps should give nil" do
    d3 = RubyLie::Algebra.new(:alg_D,3)
    d4 = RubyLie::Algebra.new(:alg_D,4)
    b3 = RubyLie::Algebra.new(:alg_B,3)
    b4 = RubyLie::Algebra.new(:alg_B,4)

    RubyLie::YoungTableau.from_highest_weight(d3.omega(2)).should be_nil
    RubyLie::YoungTableau.from_highest_weight(d3.omega(3)).should be_nil

    RubyLie::YoungTableau.from_highest_weight(d4.omega(3)).should be_nil
    RubyLie::YoungTableau.from_highest_weight(d4.omega(4)).should be_nil

    RubyLie::YoungTableau.from_highest_weight(b3.omega(3)).should be_nil

    RubyLie::YoungTableau.from_highest_weight(b4.omega(4)).should be_nil
  end
end

