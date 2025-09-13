from manim import *
config.background_color = WHITE


style = {"fill_color": BLACK, "background_stroke_color": WHITE}


class ResidualMLP(Scene):
    def construct(self):
        # Define the layers
        layer_labels = ["\\begin{bmatrix} x_1 \\\ x_2 \\end{bmatrix}", 
                        "\\begin{bmatrix} p_1 \\\ p_2 \\\ p_3 \\end{bmatrix}", 
                        "\\begin{bmatrix} h_1 \\\ h_2 \\\ h_3 \\end{bmatrix}", 
                        "\\begin{bmatrix} y_1 \\\ y_2 \\end{bmatrix}",  # LaTeX vector notation
                        "L"]
        operations = ["\\text{Linear}", "\\text{ReLU}", "\\text{Linear}", "\\text{Loss}"]
        operations2 = ["f_1", "f_2", "f_3", "f_4"]
        layer_positions = [LEFT * 6, LEFT * 3, ORIGIN, RIGHT * 3, RIGHT * 6.5] 
        vectors = []
     
        plus_symbol = MathTex("+").scale(0.9).set_color(BLACK).next_to(layer_positions[-2] + RIGHT * 0.5)
        circle_around_plus = Circle(radius=0.25, color=BLACK, fill_opacity=0.2).move_to(plus_symbol.get_center())
        res_text = MathTex("\mathbf{x}").scale(0.75).next_to(plus_symbol.get_center() + 0.05 + UP * 1, buff=0.2).set_color(BLACK)


        line_group = VGroup()
        line_group += Line(layer_positions[0] + UP * 0.75, layer_positions[0] + UP * 2, color=BLACK)
        line_group += Dot(layer_positions[0] + UP * 2, color=BLACK, radius=0.02)
        line_group += Line(layer_positions[0] + UP * 2, plus_symbol.get_center() + UP * 2, color=BLACK)
        line_group += Dot(plus_symbol.get_center() + UP * 2, color=BLACK, radius=0.02)
        line_group += Line(plus_symbol.get_center() + UP * 2, plus_symbol.get_center() + UP * circle_around_plus.get_radius(), color=BLACK)


        self.add(line_group, plus_symbol, circle_around_plus, res_text)

        # Create layer vectors
        for i, label in enumerate(layer_labels):
            vec = MathTex(label).scale(0.9).set_color(BLACK)
            vec.move_to(layer_positions[i])
            self.add(vec)
            vectors.append(vec)

        # Create operation arrows
        for i in range(len(layer_positions) - 1):
            if i == 3:
                arrow = Arrow(plus_symbol.get_center() + RIGHT * circle_around_plus.get_radius(), layer_positions[i + 1] + LEFT * 0.55, buff=0.1, stroke_width=5, color=BLACK)
            else: 
                arrow = Arrow(layer_positions[i] + RIGHT * 0.55, layer_positions[i + 1] + LEFT * 0.55, buff=0.1, stroke_width=5, color=BLACK)
            op_text = MathTex(operations[i]).scale(0.75).next_to(arrow, UP, buff=0.2).set_color(BLACK)
            un_text = MathTex(operations2[i]).scale(0.75).next_to(arrow, DOWN, buff=0.1).set_color(BLACK)
            self.add(arrow, op_text, un_text)

# To render the image, use the command in terminal:
# manim -qk MLP_diagram.py MLPDiagram --format=png
