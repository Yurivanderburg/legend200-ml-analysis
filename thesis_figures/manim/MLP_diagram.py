from manim import *
config.background_color = WHITE


style = {"fill_color": BLACK, "background_stroke_color": WHITE}


class MLPDiagram(Scene):
    def construct(self):
        # Define the layers
        layer_labels = ["\\begin{bmatrix} x_1 \\\ x_2 \\end{bmatrix}", 
                        "\\begin{bmatrix} p_1 \\\ p_2 \\\ p_3 \\end{bmatrix}", 
                        "\\begin{bmatrix} h_1 \\\ h_2 \\\ h_3 \\end{bmatrix}", 
                        "\\begin{bmatrix} y_1 \\\ y_2 \\end{bmatrix}",  # LaTeX vector notation
                        "L"]
        operations = ["\\text{Linear}", "\\text{ReLU}", "\\text{Linear}", "\\text{Loss}"]
        operations2 = ["f_1", "f_2", "f_3", "f_4"]
        layer_positions = [LEFT * 6, LEFT * 3, ORIGIN, RIGHT * 3, RIGHT * 6] 
        vectors = []

        # Create layer vectors
        for i, label in enumerate(layer_labels):
            vec = MathTex(label).scale(1.25).set_color(BLACK)
            vec.move_to(layer_positions[i])
            self.add(vec)
            vectors.append(vec)

        # Create operation arrows
        for i in range(len(layer_positions) - 1):
            arrow = Arrow(layer_positions[i] + RIGHT * 0.7, layer_positions[i + 1] + LEFT * 0.7, buff=0.1, color=BLACK)
            op_text = MathTex(operations[i]).scale(0.8).next_to(arrow, UP, buff=0.2).set_color(BLACK)
            un_text = MathTex(operations2[i]).scale(0.8).next_to(arrow, DOWN, buff=0.2).set_color(BLACK)
            self.add(arrow, op_text, un_text)

# To render the image, use the command in terminal:
# manim -qk MLP_diagram.py MLPDiagram --format=png
