module S2Geometry
  class S2Cell
    require_relative '../vendor/s2-geometry-java.jar'
    require_relative '../vendor/jsr305.jar'
    require_relative '../vendor/guava-r09.jar'

    import 'com.google.common.geometry'

    S2_LEVEL_4_X_4_METRES = 21

    def initialize(numeric_s2_id)
      @s2id = S2CellId.new(numeric_s2_id)
    end

    def self.build_from_lat_long(latitude, longitude)
      new(S2CellId.from_lat_lng(S2LatLng.from_degrees(latitude, longitude)).parent(S2_LEVEL_4_X_4_METRES).id)
    end
  end
end
