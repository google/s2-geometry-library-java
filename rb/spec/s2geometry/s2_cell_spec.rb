module S2Geometry
  require 'spec_helper'

  describe S2Cell do
    describe '.build_from_lat_long' do
      it 'should return an instance of s2cell id' do
        expect(described_class.build_from_lat_long(1, 1)).to be_an_instance_of(described_class)
      end
    end
  end
end

