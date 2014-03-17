
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
end

