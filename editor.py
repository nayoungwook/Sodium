import pygame
from atom import Atom
from openbabel import pybabel
from openbabel import pybel

wsize = (1280, 720)


cs = 50

class Button():
    def __init__(self, att, j, i):

        self.i = i
        self.j = j
        self.att = att

        self.pos = [0, 0]
        self.scale = cs
        
    def draw(self, screen, font, sc):
        j = self.j
        i = self.i
        att = self.att

        global wsize, cs

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
        global cs
        
        twm = abs(pygame.mouse.get_pos()[0] - (self.pos[0])) < 25 and\
            abs(pygame.mouse.get_pos()[1] - (self.pos[1])) < 25

        if twm:
            self.scale += (cs / 3 * 4 - self.scale) / 5

            if editor.click:
                atoms.append(Atom(self.att[self.i][self.j], editor.att))
        else:
            self.scale += (cs - self.scale) / 5
        
        
class Editor():

    def __init__(self):
        self.bgt = 0
        pygame.init()
        pygame.font.init()
        self.font = pygame.font.Font("font.ttf", 30)

        self.click = False
        self.atoms = []

        self.w = 18
        self.h = 4
        
        self.att = [
            ['H', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'He'],
            ['Li', 'Be', '', '', '', '', '', '', '', '', '', '', 'B', 'C', 'N', 'O', 'F', 'Ne'],
            ['Na', 'Mg', '', '', '', '', '', '', '', '', '', '', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
            ['K', 'Ca', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''],
        ]

        self.buttons = []
        
        for i in range(len(self.att)):
            for j in range(len(self.att[i])):
                if self.att[i][j] != '':
                    self.buttons.append(Button(self.att, j, i))

        # Atom 생성
        carbon = Atom("C", self.att)
        oxygen = Atom("O", self.att)
        hydrogen = Atom("H", self.att)
        self.atoms.append(carbon)
        self.atoms.append(oxygen)
        self.atoms.append(hydrogen)

        # Bond 생성
        carbon.add_bond(hydrogen, "-")
        carbon.add_bond(oxygen, "=")

        # Smiles 변환
        smiles_representation = self.to_smiles(carbon)
        print("Smiles 표기법:", smiles_representation)

        m = pybel.readstring("smi", smiles_representation)
        print(m)

    def to_smiles(self, atom, visited=None):
        if visited is None:
            visited = set()

        smiles = atom.at_type
        visited.add(atom)

        for neighbor, bond_type in atom.bonds:
            if neighbor not in visited:
                smiles += f"({self.to_smiles(neighbor, visited)})"
                smiles += bond_type

        return smiles

    def run(self):
        global wsize
        self.screen = pygame.display.set_mode(wsize)
        self.clock = pygame.time.Clock()
        run = True

        while run:
            self.click = False
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    run = False
                elif event.type == pygame.MOUSEBUTTONUP:
                    self.click = True

            self.tick()
            self.draw()
            
            pygame.display.flip()

            self.clock.tick(60)

        pygame.quit()

    def tick(self):
        for btn in self.buttons:
            btn.tick(self.atoms, self)
        
    def bg(self):
        self.bgt -= 1

        g = 150
        
        if self.bgt < -g:
            self.bgt = 0
        
        for i in range(15):
            for j in range(10):
                c = (200, 200, 200)
                if (i + j) % 2 == 0:
                    c = (100, 100, 100)
                    
                pygame.draw.rect(self.screen, c, (i * g + self.bgt, j * g + self.bgt, g, g))
                
    def table(self):

        sc = (wsize[0], 260)
        pygame.draw.rect(self.screen, (20, 20, 20), (0, wsize[1] - sc[1], sc[0], sc[1]))

        for btn in self.buttons:
            btn.draw(self.screen, self.font, sc)
        
    def draw(self):
        self.screen.fill((0, 55, 0))
        self.bg()
        self.table()

        for atom in self.atoms:
            atom.draw(self.screen, self.font)
