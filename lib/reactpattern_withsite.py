from reactpattern_pure import ReactPatternNew

class ReactPattern(ReactPatternNew):
    def SaveSurfaceMetalBmx(self,natom,ISbmx,FSbmx):
        for i in range(natom):
            for j in range(natom):
                if i ==1 and j == 2:

    def AddSiteMatch(self,sitematch):
        for pair in sitematch:
            for i,atominfoi in enumerate(self.atommatch):
                if atominfoi == pair[0]
