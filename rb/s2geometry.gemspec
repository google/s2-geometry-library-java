# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 's2geometry/version'

Gem::Specification.new do |s|
  s.name        = 's2geometry-jar'
  s.version     = S2Geometry::VERSION
  s.platform    = 'java'
  s.authors     = ['Sidu Ponnappa']
  s.email       = ['ckponnappa@gmail.com']
  s.homepage    = 'https://github.com/gojek-engineering/s2-geometry-library-java'
  s.summary     = %q{Google's S2 Geometry Java Library}
  s.description = %q{Gem package of Google's S2 Geometry Java Library}

  s.files         = Dir['lib/**/*.{rb,jar}']
  s.require_paths = %w(lib)
end