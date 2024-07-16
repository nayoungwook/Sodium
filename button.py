import pygame
from atom import Atom
from atom import sel_atom

class Button():
    def __init__(self, att, j, i, cs, wsize):

        self.i = i
        self.j = j
        self.att = att

        self.pos = [0, 0]
        
        self.scale = cs
        self.cs = cs
        self.wsize = wsize
        
    def draw(self, screen, font, sc):
        j = self.j
        i = self.i
        att = self.att

        cs = self.cs
        wsize = self.wsize
        
        self.pos = [wsize[0] / 2 + (cs) * (-9 + j) + cs / 2, wsize[1] - sc[1] + cs / 2 + (cs) * i + cs / 2]
        
        pygame.draw.rect(screen, (140, 140, 140), \
                         (self.pos[0] - self.scale / 2, self.pos[1] - self.scale / 2,\
                                  self.scale, self.scale))
        pygame.draw.rect(screen, (240, 240, 240), \
                             (self.pos[0] + 2 - self.scale / 2, self.pos[1] + 2 - self.scale / 2,\
                                  self.scale - 4, self.scale - 4))

        f = font.render(att[i][j], True, (0, 0, 0))
        screen.blit(f, f.get_rect(center=(wsize[0] / 2 + (cs) * (-9 + j) + cs / 2,\
                                               wsize[1] - sc[1] + cs / 2 + (cs) * i + cs / 2)))

    def tick(self, atoms, editor):
        cs = self.cs
        twm = abs(pygame.mouse.get_pos()[0] - (self.pos[0])) < 25 and\
            abs(pygame.mouse.get_pos()[1] - (self.pos[1])) < 25

        if twm:
            self.scale += (cs / 3 * 4 - self.scale) / 5

            if editor.click:
                atoms.append(Atom(self.att[self.i][self.j], editor.att))
        else:
            self.scale += (cs - self.scale) / 5
        
