import pygame
from OpenGL.GL import *

class App:
    def __init__(self):
        pygame.init()
        pygame.display.set_mode((640, 480), pygame.OPENGL | pygame.DOUBLEBUF)

        self.clock = pygame.time.Clock()

        glClearColor(0.1, 0.2, 0.2, 1)

        self.loop()

    def draw(self):
        glClear(GL_COLOR_BUFFER_BIT)
        
        pygame.display.flip()
        
    def loop(self):

        self.run = True

        while self.run:
            for event in pygame.event.get():
                if event.type == pygame.quit:
                    self.run = False


            self.draw()

            self.clock.tick(60)

        self.quit()

    def quit(self):
        pygame.quit()
