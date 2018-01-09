<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="2.18.15" minimumScale="inf" maximumScale="1e+08" hasScaleBasedVisibilityFlag="0">
  <pipe>
    <rasterrenderer opacity="1.0" alphaBand="-1" classificationMax="4" classificationMinMaxOrigin="User" band="1" classificationMin="1" type="singlebandpseudocolor">
      <rasterTransparency/>
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" clip="1">
          <item alpha="255" value="1" label="I" color="#ffffff"/>
          <item alpha="255" value="2" label="II" color="#bfccff"/>
          <item alpha="255" value="3" label="III" color="#a0e6ff"/>
          <item alpha="255" value="4" label="IV" color="#80ffff"/>
          <item alpha="255" value="5" label="V" color="#7aff93"/>
          <item alpha="255" value="6" label="VI" color="#ffff00"/>
          <item alpha="255" value="7" label="VII" color="#ffc800"/>
          <item alpha="255" value="8" label="VIII" color="#ff9100"/>
          <item alpha="255" value="9" label="IX" color="#ff0000"/>
          <item alpha="255" value="10" label="X" color="#c80000"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="-6" contrast="0"/>
    <huesaturation colorizeGreen="128" colorizeOn="0" colorizeRed="255" colorizeBlue="128" grayscaleMode="0" saturation="0" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>6</blendMode> <!-- Blending mode: multiply -->
</qgis>
