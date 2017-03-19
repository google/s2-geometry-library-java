module S2Geometry
  require 'spec_helper'
  import 'com.google.common.geometry'

  describe S2Cell do
    describe '.build_from_lat_long' do
      it 'should return an instance of s2cell id' do
        expect(described_class.build_from_lat_long(1, 1)).to be_an_instance_of(described_class)
      end

      it 'should return an instance of s2cell id which has an instance of s2idcell as attribute' do
        expect(described_class.build_from_lat_long(1, 1).s2id).to be_an_instance_of(S2CellId)
      end
    end

    describe '.convert to lat_long' do
      it 'should return an instance of latlong' do
        expect(described_class.convert_to_lat_long(3211)).to be_an_instance_of(S2LatLng)
      end
    end
  end
end
