# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'ruby_lie/version'

Gem::Specification.new do |spec|
  spec.name          = "ruby_lie"
  spec.version       = RubyLie::VERSION
  spec.authors       = ["Chris Locke"]
  spec.email         = ["project.eutopia@gmail.com"]
  spec.description   = %q{A library for dealing with Lie algebras}
  spec.summary       = %q{Handles vectors in root/weight basis, highest weight representations, and others}
  spec.homepage      = ""
  spec.license       = "MIT"

  spec.files         = `git ls-files`.split($/)
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 1.3"
  spec.add_development_dependency "rake"
  spec.add_development_dependency 'rspec'
end
