<?xml version="1.0" encoding="ISO-8859-1"?>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:text="urn:oasis:names:tc:opendocument:xmlns:text:1.0">

<xsl:template match="ul">
	<text:list>
	<xsl:apply-templates/>
	</text:list>
</xsl:template>

<xsl:template match="li">
	<text:list-item>
	<text:p><xsl:value-of select="."/></text:p>
	</text:list-item>
</xsl:template>

<xsl:template match="/">
	<xsl:apply-templates/>
</xsl:template>

</xsl:stylesheet>
