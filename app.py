import pygame
from OpenGL.GL import *
from OpenGL.GL.shaders import *
import draw

class App:
    def __init__(self):
        pygame.init()
        pygame.display.set_mode((640, 480), pygame.OPENGL | pygame.DOUBLEBUF)

        self.clock = pygame.time.Clock()

        glClearColor(0.1, 0.2, 0.2, 1)

        self.shader = self.mshader('vertex.vert', 'fragment.frag')
        self.triangle = draw.Triangle()
        
        self.loop()

    def mshader(self, vp, fp):

        vs, fs = '', ''
        
        with open(vp, 'r') as f:
            vs = f.readlines()

        with open(fp, 'r') as f:
            fs = f.readlines()

        shader = compileProgram(
            compileShader(vs, GL_VERTEX_SHADER),
            compileShader(fs, GL_FRAGMENT_SHADER)
        )

        return shader
        
    def draw(self):
        glClear(GL_COLOR_BUFFER_BIT)

        glUseProgram(self.shader)
        glBindVertexArray(self.triangle.vao)
        glDrawArrays(GL_TRIANGLES, 0, self.triangle.vertex_count)
        
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

        self.triangle.destroy()
        glDeleteProgram(self.shader)
        
        pygame.quit()
