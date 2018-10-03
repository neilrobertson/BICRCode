'''
Created on 27 Nov 2014

@author: johncole
'''

from bed.treatment import Bed as BedFile
import math

def getFractionalTagsInRegion(chrm, start, end, tagsInRegion):
    values = []
    for tag in tagsInRegion:
        assert tag.end>tag.start       
        containedlength = min(tag.end,end) - max(tag.start,start)
        length = tag.end - tag.start      
        values.append(float(containedlength) / float(length)) 
    return values
        

class ChIP():
    
    #Constructor
    def __init__(self,filename,extend=None):
        self.data = BedFile(filename,extend=extend)
        self.getTagsInRegion = self.data.getIntervalsInRange
        
    #Gets the reads per million for the library:
    def getLibraryNormalizationFactor(self):
        return self.data.count / 1000000

    #Gets normalised tag counts for a region - using a variety of normalisation strategies
    def normaliseTags(self,chrm, start, end, tagsInRegion):
        
        assert end > start, str(end) + " > " + str(start)
        fractionalTagsInRegion = getFractionalTagsInRegion(chrm, start, end, tagsInRegion)
        
        Tags = math.fsum(fractionalTagsInRegion)
        TagsRegionNormalised = Tags / float(end-start)
        TagsLibraryNormalised = Tags/ self.getLibraryNormalizationFactor()
        TagsRegionLibraryNormalised = TagsRegionNormalised / self.getLibraryNormalizationFactor()
        
        return [Tags,TagsRegionNormalised,TagsLibraryNormalised,TagsRegionLibraryNormalised]

class ChIPWithControl():
    
    #Constructor
    def __init__(self,ChIPFilename, controlFilename, extend=None):
        self.ChIPData = BedFile(ChIPFilename,extend=extend)
        self.ControlData = BedFile(controlFilename,extend=extend)
    
    #Gets the tag counts for both ChIP and control datasets
    def getTagsInRegion(self,chrm, start, stop):
        # store a tuple of sample/cntrl values
        return self.ChIPData.getIntervalsInRange(chrm,start,stop),self.ControlData.getIntervalsInRange(chrm, start, stop) 
     
    #Gets the reads per million for the library:
    def getLibraryNormalizationFactor(self,data):
        return data.count / 1000000
    
    #Gets normalised tag counts for a region - using a variety of normalisation strategies
    def normaliseTags(self,chrm, start, end, chipTagsInRegion, controlTagsInRegion):
        
        assert end > start, str(end) + " > " + str(start)
        chipFractionalTagsInRegion = getFractionalTagsInRegion(chrm, start, end, chipTagsInRegion)
        controlFractionalTagsInRegion = getFractionalTagsInRegion(chrm, start, end, controlTagsInRegion)


        chipTags = math.fsum(chipFractionalTagsInRegion)
        chipTagsRegionNormalised = chipTags / float(end-start)
        chipTagsLibraryNormalised = chipTags/ self.getLibraryNormalizationFactor(self.ChIPData)
        chipTagsRegionLibraryNormalised = chipTagsRegionNormalised / self.getLibraryNormalizationFactor(self.ChIPData)
        
        controlTags = math.fsum(controlFractionalTagsInRegion)
        controlTagsRegionNormalised = controlTags / float(end-start)
        controlTagsLibraryNormalised = controlTags/ self.getLibraryNormalizationFactor(self.ControlData)
        controlTagsRegionLibraryNormalised = controlTagsRegionNormalised / self.getLibraryNormalizationFactor(self.ControlData)
        
        diffTags = chipTags - controlTags
        diffTagsRegionNormalised = chipTagsRegionNormalised - controlTagsRegionNormalised
        diffTagsLibraryNormalised = chipTagsLibraryNormalised - controlTagsLibraryNormalised
        diffTagsRegionLibraryNormalised = chipTagsRegionLibraryNormalised - controlTagsRegionLibraryNormalised
        
        if (controlTags > 0):
            ratioTags = chipTags / controlTags
            ratioTagsRegionNormalised = chipTagsRegionNormalised / controlTagsRegionNormalised
            ratioTagsLibraryNormalised = chipTagsLibraryNormalised / controlTagsLibraryNormalised
            ratioTagsRegionLibraryNormalised = chipTagsRegionLibraryNormalised / controlTagsRegionLibraryNormalised
        else:
            ratioTags = 0.0
            ratioTagsRegionNormalised = 0.0
            ratioTagsLibraryNormalised = 0.0
            ratioTagsRegionLibraryNormalised = 0.0

        
        return chipTags,chipTagsRegionNormalised,chipTagsLibraryNormalised,chipTagsRegionLibraryNormalised,controlTags,controlTagsRegionNormalised,controlTagsLibraryNormalised,controlTagsRegionLibraryNormalised,diffTags,diffTagsRegionNormalised,diffTagsLibraryNormalised,diffTagsRegionLibraryNormalised,ratioTags,ratioTagsRegionNormalised,ratioTagsLibraryNormalised,ratioTagsRegionLibraryNormalised
    
    
    