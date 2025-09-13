from manim import *


#style = {"fill_color": BLACK, "background_stroke_color": WHITE}
config.background_color = WHITE

class AttentionPattern(Scene):
    def construct(self):

        m1 = Matrix([["q_1 k_1", "q_2 k_1", "q_3 k_1"],
                     ["q_1 k_2", "q_2 k_2", "q_3 k_2"], 
                     ["q_1 k_3", "q_2 k_3", "q_3 k_3"]]).set_color(BLACK)
        m1.add(SurroundingRectangle(m1.get_columns()[1]).set_color(GREEN))
        #m1.rotate(+30 * DEGREES, axis=UP)

        # Create label below column 2
        column_label = MathTex("\\sum_{i=1}^3 \\frac{q_2 k_i}{\\sqrt{d_k}} = 1 \\implies \\mathrm{softmax}")  # Add your LaTeX text here
        column_label.scale(0.7)  # Make it a bit smaller
        column_label.set_color(BLACK)
        column_label.next_to(m1.get_columns()[1], DOWN * 1.2)  # Position it below the last cell

        self.add(m1, column_label)
