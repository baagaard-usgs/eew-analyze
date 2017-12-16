<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="2.18.15" minimumScale="inf" maximumScale="1e+08" hasScaleBasedVisibilityFlag="0">
  <pipe>
    <rasterrenderer opacity="0.50" alphaBand="-1" classificationMax="1000" classificationMinMaxOrigin="User" band="1" classificationMin="0" type="singlebandpseudocolor">
      <rasterTransparency/>
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" clip="0">
          <item alpha="255" value="0" label="0" color="#fafafa"/>
          <item alpha="255" value="0.1" label="0.1" color="#c9c9c9"/>
          <item alpha="255" value="1" label="1" color="#989898"/>
          <item alpha="255" value="10" label="10" color="#676767"/>
          <item alpha="255" value="100" label="100" color="#363636"/>
          <item alpha="255" value="1000" label="1e+03" color="#050505"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="-6" contrast="0"/>
    <huesaturation colorizeGreen="128" colorizeOn="0" colorizeRed="255" colorizeBlue="128" grayscaleMode="0" saturation="0" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>6</blendMode> <!-- Blending mode: multiply -->
</qgis>
