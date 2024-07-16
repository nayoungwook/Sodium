import pygame
import math

sel_atom = None
bond_atom = None

class Vector:

    def __init__(self, x, y):
        self.x = x
        self.y = y


class Atom:

    def __init__(self, _type, att):
        self.at_type = _type
        self.att = att
        self.position = Vector(1280 / 2, 720 / 2)
        self.bonds = []
        self.vis = False

    def draw_bond(self, screen):
        for bond in self.bonds:
#            print(bond[1])
            if bond[1] == '-':
                pygame.draw.line(screen, (255, 255, 240), (bond[0].position.x, bond[0].position.y), (self.position.x, self.position.y), 5)
            elif bond[1] == '=':
                pygame.draw.line(screen, (255, 255, 240), (bond[0].position.x - 5, bond[0].position.y - 5), (self.position.x - 5, self.position.y - 5), 3)
                pygame.draw.line(screen, (255, 255, 240), (bond[0].position.x + 5, bond[0].position.y + 5), (self.position.x + 5, self.position.y + 5), 3)
                
    def draw(self, screen, font):
        pygame.draw.circle(screen, (150, 255, 155), (self.position.x, self.position.y), 20)
        pygame.draw.circle(screen, (100, 200, 55), (self.position.x, self.position.y), 20, 6)
        f = font.render(self.at_type, True, (0, 0, 0))
        screen.blit(f, f.get_rect(center=(self.position.x, self.position.y)))


    def add_bond(self, atom, bond_type):
        if isinstance(atom, Atom) and isinstance(bond_type, str):
            self.bonds.append((atom, bond_type))
            atom.bonds.append((self, bond_type))
        else:
            raise ValueError("Invalid Atom or bond type")

    def __repr__(self):
        return f"Atom({self.at_type})"

    def tick(self, atoms):
        global sel_atom, bond_atom
        
        if math.sqrt( (self.position.x - pygame.mouse.get_pos()[0]) ** 2 + (self.position.y - pygame.mouse.get_pos()[1]) **2) <= 30:
            if pygame.mouse.get_pressed()[0]:
                if sel_atom == None:
                    sel_atom = self
                
            if pygame.mouse.get_pressed()[2]:
                if bond_atom == None:
                    bond_atom = self
                
        if not pygame.mouse.get_pressed()[0]:
            if sel_atom == self:
                sel_atom = None
                    
        if not pygame.mouse.get_pressed()[2]:
            if bond_atom == self:
                for atom in atoms:
                    if atom != self:
                        if math.sqrt( (atom.position.x - pygame.mouse.get_pos()[0]) ** 2 + (atom.position.y - pygame.mouse.get_pos()[1]) ** 2 ) <= 30:
                            self.add_bond(atom, '-')
                
                bond_atom = None
                    
        if sel_atom == self:
            self.position.x = pygame.mouse.get_pos()[0]
            self.position.y = pygame.mouse.get_pos()[1]
                    
        for atom in atoms:
            if atom != self:
                if math.sqrt( (self.position.x - atom.position.x) ** 2 + (self.position.y - atom.position.y) ** 2) < 30:
                    if self.position.x == atom.position.x and self.position.y == atom.position.y:
                        self.position.x += 1
                        atom.position.x -= 1
                        self.position.y -= 1
                        atom.position.y += 1
                    
                    self.position.x += (self.position.x - atom.position.x) / 5
                    self.position.y += (self.position.y - atom.position.y) / 5
        
class Bond:
    def __init__(self, atom1, atom2, bond_type):
        if isinstance(atom1, Atom) and isinstance(atom2, Atom) and isinstance(bond_type, str):
            self.atom1 = atom1
            self.atom2 = atom2
            self.bond_type = bond_type
        else:
            raise ValueError("Invalid Atom or bond type")

    def __repr__(self):
        return f"Bond({self.atom1.symbol}-{self.bond_type}-{self.atom2.symbol})"
