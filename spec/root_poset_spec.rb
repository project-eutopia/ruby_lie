require 'spec_helper'

describe RubyLie::RootPoset do
  ALGEBRAS.each do |algebra|
    context "#{algebra.to_latex}" do

      it "highest root from poset agrees with -alpha_0" do
        poset = algebra.root_poset
        expect(poset.highest_root).to be == -algebra.alpha(0)
      end

    end
  end
end
