import pygame
from atom import Atom
from atom import sel_atom
from button import Button
from openbabel import pybel

wsize = (1280, 720)

cs = 50
        
class Editor():

    def __init__(self):
        self.bgt = 0
        pygame.init()
        pygame.font.init()
        self.font = pygame.font.Font("font.ttf", 30)

        self.click = False
        self.atoms = []

        self.smi = ''
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
                    self.buttons.append(Button(self.att, j, i, cs, wsize))

    def cr_m(self): # create model
        
        for atom in self.atoms:
            atom.vis = False
            
        for atom in self.atoms:
            if atom.vis:
                continue

            self.smi = self.to_smiles(atom)
            print('created : ', self.smi)

            if self.smi != None:
                mol = pybel.readstring('smi', self.smi)
                mol.make3D()
                print(mol.write('sdf'))
                mol.draw()

        
    def to_smiles(self, atom, visited=None):
        if visited is None:
            visited = set()

        if atom.at_type == 'H':
            return
            
        smiles = atom.at_type
        visited.add(atom)
        atom.vis = True

        for neighbor, bond_type in atom.bonds:
            if neighbor not in visited:
                if neighbor.at_type != 'H':
                    if bond_type == '=':
                        smiles += '1'

                    smiles += f"{self.to_smiles(neighbor, visited)}"

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
                elif event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_SPACE:
                        self.cr_m()
                        
            self.tick()
            self.draw()
            
            pygame.display.flip()

            self.clock.tick(60)

        pygame.quit()

    def tick(self):
        for btn in self.buttons:
            btn.tick(self.atoms, self)

        for atom in self.atoms:
            atom.tick(self.atoms)
            atom.vis = False
        
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

        from atom import bond_atom
        
        if bond_atom != None:
            pygame.draw.line(self.screen, (255, 255, 240), (bond_atom.position.x, bond_atom.position.y), pygame.mouse.get_pos(), 5)

        for atom in self.atoms:
            atom.draw_bond(self.screen)
            
        for atom in self.atoms:
            atom.draw(self.screen, self.font)
