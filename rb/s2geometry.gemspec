# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 's2geometry/version'

Gem::Specification.new do |spec|
  spec.name        = 's2geometry-jar'
  spec.version     = S2Geometry::VERSION
  spec.platform    = 'java'
  spec.authors     = ['Sidu Ponnappa']
  spec.email       = ['ckponnappa@gmail.com']
  spec.homepage    = 'https://github.com/gojek-engineering/s2-geometry-library-java'
  spec.summary     = %q{Google's S2 Geometry Java Library}
  spec.description = %q{Gem package of Google's S2 Geometry Java Library}

  spec.files         = Dir['lib/**/*.{rb,jar}']
  spec.require_paths = %w(lib)
  spec.add_dependency "jbundler"
  spec.add_development_dependency "simplecov", "~> 0.11"
  spec.add_development_dependency "bundler", "~> 1.8"
  spec.add_development_dependency "rake"
  spec.add_development_dependency "rspec", "~> 3"
  spec.add_development_dependency 'ruby-maven', '~> 3.3'

  spec.add_dependency 'sentry-raven', '~> 1.2.3' 
end
