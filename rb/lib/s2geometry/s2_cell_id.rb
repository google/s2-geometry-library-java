module S2Geometry
  class S2CellId
    import 'com.google.common.geometry'

    S2_LEVEL_4_X_4_METRES = 21

    def initialise(numeric_s2_id)
      @s2id = S2CellId.new(numeric_s2id)
    end

    def self.build_from_lat_long(latitude, longitude)
      new(S2CellId.from_lat_lng(S2LatLng.from_degrees(latitude, longitude)).parent(S2_LEVEL_4_X_4_METRES).id)
    end
  end
end

