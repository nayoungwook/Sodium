import pygame
import math

class Atom:

    def __init__(self, _type, att):
        self.at_type = _type
        self.att = att
        self.position = [1280 / 2, 720 / 2]
        self.bonds = []

    def draw(self, screen, font):
        pygame.draw.circle(screen, (150, 100, 255), self.position, 30)
        pygame.draw.circle(screen, (100, 60, 155), self.position, 30, 6)
        f = font.render(self.at_type, True, (0, 0, 0))
        screen.blit(f, f.get_rect(center=(self.position[0], self.position[1])))

    def add_bond(self, atom, bond_type):
        if isinstance(atom, Atom) and isinstance(bond_type, str):
            self.bonds.append((atom, bond_type))
            atom.bonds.append((self, bond_type))
        else:
            raise ValueError("Invalid Atom or bond type")

    def __repr__(self):
        return f"Atom({self.symbol})"


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
