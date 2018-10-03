<xsl:template name="strip-tags">
<xsl:param name="text"/>
<xsl:choose>
<xsl:when test="contains($text, '&lt;')">
<xsl:value-of select="substring-before($text, '&lt;')"/>
<xsl:call-template name="strip-tags">
<xsl:with-param name="text" select="substring-after($text, '&gt;')"/>
</xsl:call-template>
</xsl:when>
<xsl:otherwise>
<xsl:value-of select="$text"/>
</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:call-template name="strip-tags">
<xsl:with-param name="text" select="Profile"/>
</xsl:call-template>
