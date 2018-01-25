<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="2.18.15" minimumScale="inf" maximumScale="1e+08" hasScaleBasedVisibilityFlag="0">
  <pipe>
    <rasterrenderer opacity="1" alphaBand="-1" classificationMax="1.5" classificationMinMaxOrigin="User" band="1" classificationMin="-1.5" type="singlebandpseudocolor">
      <rasterTransparency/>
      <rastershader>
        <colorrampshader colorRampType="DISCRETE" clip="0">
          <item alpha="255" value="0.5" label="True Negative" color="#abd9e9"/>
          <item alpha="255" value="1.5" label="False Negative" color="#d7191c"/>
          <item alpha="255" value="2.5" label="False Positive" color="#fdae61"/>
          <item alpha="255" value="3.5" label="True Positive" color="#2c7bb6"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0"/>
    <huesaturation colorizeGreen="128" colorizeOn="0" colorizeRed="255" colorizeBlue="128" grayscaleMode="0" saturation="0" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>6</blendMode> <!-- Blending mode: multiply -->
</qgis>
