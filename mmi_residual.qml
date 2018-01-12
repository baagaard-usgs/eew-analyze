<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="2.18.14" minimumScale="inf" maximumScale="1e+08" hasScaleBasedVisibilityFlag="0">
  <pipe>
    <rasterrenderer opacity="1" alphaBand="0" classificationMax="1.5" classificationMinMaxOrigin="User" band="1" classificationMin="-1.5" type="singlebandpseudocolor">
      <rasterTransparency/>
      <rastershader>
        <colorrampshader colorRampType="DISCRETE" clip="0">
          <item alpha="255" value="-999" label="No Data" color="#999999"/>
          <item alpha="255" value="-1.5" label="&lt;= -1.5" color="#0571b0"/>
          <item alpha="255" value="-1.25" label="-1.5 - -1.25" color="#348dc0"/>
          <item alpha="255" value="-1" label="-1.25 - -1.0" color="#63a9cf"/>
          <item alpha="255" value="-0.75" label="-1.0 - -0.75" color="#92c5de"/>
          <item alpha="255" value="-0.5" label="-0.75 - -0.50" color="#b4d6e7"/>
          <item alpha="255" value="-0.25" label="-0.50 - -0.25" color="#d6e7ef"/>
          <item alpha="255" value="0.25" label="-0.25 - 0.25" color="#f7f7f7"/>
          <item alpha="255" value="0.5" label="0.25 - 0.50" color="#f6dcd0"/>
          <item alpha="255" value="0.75" label="0.50 - 0.75" color="#f5c1a9"/>
          <item alpha="255" value="1" label="0.75 - 1.0" color="#f4a582"/>
          <item alpha="255" value="1.25" label="1.0 - 1.25" color="#e66e61"/>
          <item alpha="255" value="1.5" label="1.25 - 1.50" color="#d83740"/>
          <item alpha="255" value="inf" label="> 1.50" color="#ca0020"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0"/>
    <huesaturation colorizeGreen="128" colorizeOn="0" colorizeRed="255" colorizeBlue="128" grayscaleMode="0" saturation="0" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>6</blendMode> <!-- Blending mode: multiply -->
</qgis>
