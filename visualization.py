from manim import *
import math

class Chapter0PrologueScene(Scene):
    def construct(self):
        title_Page = Text("Euler's Sum of Powers Conjecture", font_size=48)
        title_Page_subtitle = Paragraph("Exploring Euler’s Sum of Powers Through Exact Computation\nModular Arithmetic, Sumsets, and Search Design", 
                    font_size=25, alignment="center").next_to(title_Page, DOWN, buff=0.7)
        combined_object = VGroup(title_Page, title_Page_subtitle).center()
        self.play(Write(combined_object))
        self.wait(3)
        self.play(FadeOut(combined_object))

        title = Text("Chapter 0: Prologue", font_size=46).to_edge(UP)
        self.play(Write(title), run_time=1.4)
        self.wait(5)

        euler = Text("Leonhard Euler (1769)", font_size=32, color=BLUE_C)
        claim = Text(
            "Conjecture: need at least k terms",
            font_size=30,
            color=WHITE,
        ).next_to(euler, DOWN, buff=0.35)
        statement = MathTex(
            r"a_1^k + a_2^k + \cdots + a_n^k = b^k,\quad k>2",
            font_size=30,
            color=YELLOW,
        ).next_to(claim, DOWN, buff=0.35)

        block1 = VGroup(euler, claim, statement).move_to(UP * 0.7)

        self.play(Write(euler), run_time=0.9)
        self.play(FadeIn(claim, shift=UP * 0.2), run_time=0.8)
        self.play(Write(statement), run_time=1.4)
        self.wait(5)

        busted = Text("1966: counterexample found", font_size=30, color=RED_C).move_to(DOWN * 1.25)
        counterexample = MathTex(
            r"27^5 + 84^5 + 110^5 + 133^5 = 144^5",
            font_size=34,
            color=GREEN_C,
        ).next_to(busted, DOWN, buff=0.35)
        note = Text("Only 4 terms on the left for k=5", font_size=26, color=GRAY_A).next_to(
            counterexample, DOWN, buff=0.35
        )
        block2 = VGroup(busted, counterexample, note)

        #self.play(FadeOut(block1, shift=UP * 0.3), run_time=0.7)
        self.play(FadeIn(block2, shift=DOWN * 0.2), run_time=1.1)
        self.wait(5)

        def nondec_count(n: int, r: int) -> int:
            return math.comb(n + r - 1, r)

        max_b = 250
        lhs_r = 4
        rhs_8 = 4
        rhs_9 = 5
        lhs_table = nondec_count(max_b, lhs_r)
        rhs_per_b_m8 = nondec_count(max_b, rhs_8)
        rhs_per_b_m9 = nondec_count(max_b, rhs_9)
        rhs_total_m8 = sum(nondec_count(b, rhs_8) for b in range(2, max_b + 1))
        rhs_total_m9 = sum(nondec_count(b, rhs_9) for b in range(2, max_b + 1))

        q = Text("What happens at k = 9?", font_size=44, color=YELLOW).move_to(UP * 2)
        search = Text("MITM split from code: r=floor(m/2), tuples are nondecreasing", font_size=24, color=GRAY_A)
        search.next_to(q, DOWN, buff=0.25)
        formula = MathTex("Estimate Tuple Count (N,r) = C(N+r-1, r)", font_size=26, color=TEAL_B)
        formula.next_to(search, DOWN, buff=0.3)

        lhs_box = RoundedRectangle(width=5.2, height=1.9, corner_radius=0.16, color=BLUE_C)
        lhs_box.shift(LEFT * 3.4 + DOWN * 0.6)
        lhs_title = Text("Precompute LHS table (r=4)", font_size=24, color=BLUE_A).move_to(lhs_box.get_top() + DOWN * 0.3)
        lhs_count = Text(f"C(253,4) = {lhs_table:,}", font_size=28, color=WHITE).move_to(lhs_box.get_center() + DOWN * 0.15)
        lhs_group = VGroup(lhs_box, lhs_title, lhs_count)

        rhs_box = RoundedRectangle(width=5.2, height=1.9, corner_radius=0.16, color=YELLOW_D)
        rhs_box.shift(RIGHT * 3.4 + DOWN * 0.6)
        rhs_title = Text("Scan RHS per b (m-r terms)", font_size=24, color=YELLOW_A).move_to(rhs_box.get_top() + DOWN * 0.3)
        rhs_count_m8 = Text(f"m=8: C(253,4) = {rhs_per_b_m8:,}", font_size=24, color=WHITE).move_to(
            rhs_box.get_center() + UP * 0.12
        )
        rhs_count_m9 = Text(f"m=9: C(254,5) = {rhs_per_b_m9:,}", font_size=24, color=RED_C).move_to(
            rhs_box.get_center() + DOWN * 0.32
        )
        rhs_group = VGroup(rhs_box, rhs_title, rhs_count_m8, rhs_count_m9)

        bridge = Arrow(lhs_group.get_right(), rhs_group.get_left(), buff=0.2, color=GRAY_B, stroke_width=5)
        ratio_line = Text(f"Per-b RHS growth at b=250: x{rhs_per_b_m9 / rhs_per_b_m8:.1f}", font_size=28, color=RED_B)
        ratio_line.next_to(VGroup(lhs_group, rhs_group), DOWN * 0.5, buff=0.25)

        totals = Text(
            f"Total RHS tuples for b=2...250:  m=8 -> {rhs_total_m8:,}    m=9 -> {rhs_total_m9:,}",
            font_size=22,
            color=GRAY_A,
        ).to_edge(DOWN)

        explosion = VGroup(q, search, formula, lhs_group, rhs_group, bridge, ratio_line, totals)
        self.play(FadeOut(block1, shift=UP * 0.3), run_time=0.7)
        self.play(FadeOut(block2, shift=UP * 0.25), run_time=0.7)
        self.play(Write(q), run_time=0.7)
        self.wait(3)
        self.play(FadeIn(search, shift=UP * 0.2), FadeIn(formula, shift=UP * 0.2), run_time=0.8)
        self.wait(5)
        self.play(FadeIn(lhs_group, shift=RIGHT * 0.2), run_time=0.8)
        self.play(FadeIn(rhs_group, shift=LEFT * 0.2), GrowArrow(bridge), run_time=1.0)
        self.play(FadeIn(ratio_line, shift=UP * 0.15), run_time=0.6)
        self.play(FadeIn(totals, shift=UP * 0.1), run_time=0.6)
        self.wait(5)

        left = RoundedRectangle(width=4.6, height=2.0, corner_radius=0.2, color=TEAL_C)
        left_title = Text("Exact Research Solver", font_size=28, color=TEAL_A).move_to(
            left.get_top() + DOWN * 0.25
        )
        left_text = Text(
            "Memory-bounded\nwithout losing correctness",
            font_size=24,
            line_spacing=0.8,
        ).move_to(left).shift(DOWN*0.025)
        left_group = VGroup(left, left_title, left_text).shift(LEFT * 3.2 + DOWN * 0.25)

        right = RoundedRectangle(width=4.6, height=2.0, corner_radius=0.2, color=PURPLE_C)
        right_title = Text("AI Systems Instincts", font_size=28, color=PURPLE_A).move_to(
            right.get_top() + DOWN * 0.25
        )
        right_text = Text(
            "Batch Size vs Memory\nIndexing vs Streaming\nExact vs Approximate",
            font_size=24,
            line_spacing=0.8,
        ).move_to(right).shift(DOWN*0.15)
        right_group = VGroup(right, right_title, right_text).shift(RIGHT * 3.2 + DOWN * 0.25)

        arrow = Arrow(
            left_group.get_right() + RIGHT * 0.1,
            right_group.get_left() + LEFT * 0.1,
            buff=0.1,
            color=GRAY_B,
            stroke_width=5,
        )

        self.play(FadeOut(explosion, shift=UP * 0.2), run_time=0.7)
        self.play(FadeIn(left_group, shift=RIGHT * 0.2), run_time=0.9)
        self.play(FadeIn(right_group, shift=LEFT * 0.2), run_time=0.9)
        self.play(GrowArrow(arrow), run_time=0.7)

        outro = Text("This goes beyond number theory into systems design", font_size=28, color=YELLOW_B)
        outro.to_edge(DOWN)
        self.play(FadeIn(outro, shift=UP * 0.15), run_time=0.8)
        self.wait(5)

class Chapter1ProblemScene_BruteForce(Scene):
    def construct(self):
        title_Page = Text("Chapter 1: Problem", font_size=48)
        title_Page_subtitle = Paragraph("Brute Force Computation and Meet in the Middle Optimization", 
                    font_size=25, alignment="center").next_to(title_Page, DOWN, buff=0.7)
        combined_object = VGroup(title_Page, title_Page_subtitle).center()
        self.play(Write(combined_object))
        self.wait(3)
        self.play(FadeOut(combined_object))

        title = Text("Chapter 1: Problem Setup", font_size=44).to_edge(UP)
        self.play(Write(title), run_time=1.2)

        problem = MathTex(
            r"a_1^5 + a_2^5 + \cdots + a_m^5 = b^5,\quad 1 \le a_i \le b \le N",
            font_size=34,
            color=YELLOW,
        ).next_to(title, DOWN, buff=0.45)
        objective = Text(
            "Goal: exact search for integer solutions under a finite bound N",
            font_size=26,
            color=GRAY_A,
        ).next_to(problem, DOWN, buff=0.3)
        self.play(Write(problem), FadeIn(objective, shift=UP * 0.15), run_time=1.0)
        self.wait(7)
        known_identity = MathTex(
            r"27^5 + 84^5 + 110^5 + 133^5 = 144^5",
            font_size=34,
            color=YELLOW,
        ).next_to(objective, DOWN, buff=0.35)
        self.play(Write(known_identity))
        self.wait(10)

        brute_box = RoundedRectangle(width=6, height=1.8, corner_radius=0.15, color=RED_C)
        brute_title = Text("Brute Force", font_size=30, color=RED_A).move_to(
            brute_box.get_top() + DOWN * 0.3
        )
        brute_text = Paragraph(
            "Enumerate all tuples\nand test equality directly",
            font_size=24,
            alignment="center",
            line_spacing=0.8,
        ).move_to(brute_box.get_center() + DOWN * 0.1)
        brute_formula = MathTex(r"\text{Work} \sim O(N^m)", font_size=30, color=RED_B).next_to(
            brute_box, DOWN, buff=0.2
        )
        brute_group = VGroup(brute_box, brute_title, brute_text, brute_formula
                             ).next_to(known_identity, DOWN * 1.2)
        self.play(FadeIn(brute_group), run_time=0.9)
        self.wait(10)

        N = 250
        m = 9
        r = 4
        s = m - r
        brute_nondec = math.comb(N + m - 1, m)
        lhs_count = math.comb(N + r - 1, r)
        rhs_per_b = math.comb(N + s - 1, s)
        rhs_total = sum(math.comb(b + s - 1, s) for b in range(2, N + 1))
        mitm_units = lhs_count + rhs_total
        speedup = brute_nondec / mitm_units

        numbers_title = Text("Concrete Scale Example (N=250, m=9, r=4)", font_size=32, color=YELLOW_B).next_to(brute_group, DOWN * 1.7)
        line1 = Text(f"Brute nondecreasing tuples: C(258,9) = {brute_nondec:,}", font_size=26, color=RED_C).next_to(numbers_title, DOWN * 1.2)
        self.play(FadeIn(numbers_title, shift=DOWN * 0.2), run_time=0.7)
        self.play(Write(line1), run_time=0.8)
        self.wait(10)

class Chapter1ProblemScene_MITM(Scene):
    def construct(self):
        self.wait(3.5)
        title = Text("Meet-In-The-Middle Solving", font_size=44).to_edge(UP)
        self.play(Write(title), run_time=1.0)
        self.wait(3)

        intro = Text(
            "Split the 8-term search into two 4-term halves",
            font_size=28,
            color=YELLOW_B,
        ).next_to(title, DOWN, buff=0.35)
        self.play(FadeIn(intro, shift=UP * 0.15), run_time=0.7)

        lhs_box = RoundedRectangle(width=5.2, height=2.2, corner_radius=0.15, color=TEAL_C)
        lhs_box.shift(LEFT * 3.3)
        lhs_title = Text("Precompute Index", font_size=28, color=TEAL_A).move_to(
            lhs_box.get_top() + DOWN * 0.3
        )
        lhs_text = MathTex(
            r"\text{All sums of} \\ x_1^5 + x_2^5 + x_3^5 + x_4^5... \\ \text{Store witness tuples compactly}",
            font_size=30,
            tex_environment="gather*"
        ).move_to(lhs_box.get_center() + DOWN * 0.1)
        lhs_group = VGroup(lhs_box, lhs_title, lhs_text)

        rhs_box = RoundedRectangle(width=5.2, height=2.2, corner_radius=0.15, color=BLUE_C)
        rhs_box.shift(RIGHT * 3.3)
        rhs_title = Text("Lookup Phase", font_size=28, color=BLUE_A).move_to(
            rhs_box.get_top() + DOWN * 0.3
        )
        rhs_text = MathTex(
            r"\text{For each } b^k, \\ \text{scan complementary 4-tuples}  \\ \text{and query the index}",
            font_size=30,
            tex_environment="gather*"
        ).move_to(rhs_box.get_center() + DOWN * 0.1)
        rhs_group = VGroup(rhs_box, rhs_title, rhs_text)

        flow = Arrow(lhs_group.get_right(), rhs_group.get_left(), buff=0.2, color=GRAY_B, stroke_width=5)
        db_note = Text(
            "Database analogy: build an index to later complete fast lookups",
            font_size=24,
            color=GRAY_A,
        ).next_to(VGroup(lhs_group, rhs_group), DOWN, buff=0.35)

        self.play(FadeIn(lhs_group, shift=RIGHT * 0.2), run_time=0.8)
        self.wait(4)
        self.play(FadeIn(rhs_group, shift=LEFT * 0.2), GrowArrow(flow), run_time=0.8)
        self.wait(4)
        self.play(FadeIn(db_note, shift=UP * 0.1), run_time=0.6)
        self.wait(4)

        eng_title = Text("Engineering Improvements", font_size=34, color=YELLOW_B).to_edge(UP)
        chip1 = RoundedRectangle(width=4.5, height=1.2, corner_radius=0.12, color=GREEN_C).shift(UP * 1.5)
        chip1_t = Text("Packed witness tuples", font_size=24).move_to(chip1)
        chip2 = RoundedRectangle(width=4.5, height=1.2, corner_radius=0.12, color=GREEN_C)
        chip2_t = Text("Modular residue filters", font_size=24).move_to(chip2)
        chip3 = RoundedRectangle(width=4.5, height=1.2, corner_radius=0.12, color=GREEN_C).shift(DOWN * 1.5)
        chip3_t = Text("Skip impossible candidates", font_size=24).move_to(chip3)

        self.play(
            FadeOut(VGroup(lhs_group, rhs_group, flow, db_note), shift=UP * 0.2),
            Transform(title, eng_title),
            FadeOut(intro),
            run_time=0.8,
        )

        self.play(FadeIn(chip1, shift=UP * 0.1), Write(chip1_t), run_time=0.6)
        self.play(FadeIn(chip2, shift=UP * 0.1), Write(chip2_t), run_time=0.6)
        self.play(FadeIn(chip3, shift=UP * 0.1), Write(chip3_t), run_time=0.6)
        self.wait(3)

        wall_title = Text("k = 8 Memory Wall", font_size=40, color=RED_C).to_edge(UP)
        formula = Text(
            "Nondecreasing r-tuples: C(N+r-1, r)  ~  O(N^r)",
            font_size=30,
            color=YELLOW_B,
        ).next_to(wall_title, DOWN, buff=0.25)
        split_note = Text("For a 4-tuple split: table size grows like N^4", font_size=28, color=RED_A).next_to(
            formula, DOWN, buff=0.25
        )

        bar1 = Rectangle(width=1.6, height=1, color=BLUE_C, fill_opacity=0.6).shift(LEFT * 2.2 + DOWN * 1.5)
        bar2 = Rectangle(width=1.6, height=4, color=RED_C, fill_opacity=0.6).shift(RIGHT * 2.2 + DOWN * 0.0)
        b1_label = Text("N", font_size=28).next_to(bar1, DOWN, buff=0.2)
        b2_label = Text("2N", font_size=28).next_to(bar2, DOWN, buff=0.2)
        growth = Text("Doubling N -> about 16x memory", font_size=30, color=RED_B).to_edge(DOWN)

        bar_graphs = VGroup(bar1, bar2, b1_label, b2_label).next_to(split_note, DOWN * 0.5, buff=0.25)

        self.play(
            FadeOut(VGroup(chip1, chip1_t, chip2, chip2_t, chip3, chip3_t), shift=UP),
            Transform(title, wall_title),
            run_time=0.8,
        )
        self.play(FadeIn(formula, shift=UP * 0.12), FadeIn(split_note, shift=UP * 0.12), run_time=0.7)
        self.play(FadeIn(bar_graphs))
        self.play(FadeIn(growth, shift=UP * 0.15), run_time=0.7)
        self.wait(4)

        pivot_title = Text("From Math Problem to Systems Problem", font_size=36, color=YELLOW_B).to_edge(UP)
        q = Text(
            "How do we keep search exact with bounded peak memory?",
            font_size=30,
            color=WHITE,
        ).next_to(pivot_title, DOWN, buff=0.4)

        t1 = RoundedRectangle(width=3.9, height=1.25, corner_radius=0.12, color=TEAL_C).shift(UP * 1.20)
        t1_t = Text("Precompute vs Streaming", font_size=23).move_to(t1)
        t2 = RoundedRectangle(width=3.9, height=1.25, corner_radius=0.12, color=TEAL_C).shift(DOWN * 0.4)
        t2_t = Text("Batch Size vs Memory", font_size=23).move_to(t2)
        t3 = RoundedRectangle(width=3.9, height=1.25, corner_radius=0.12, color=TEAL_C).shift(DOWN * 2)
        t3_t = Text("Exactness vs Approx.", font_size=23).move_to(t3)
        close = Text(
            "Designing around memory bottlenecks enables efficient exact computation.",
            font_size=24,
            color=GRAY_A,
        ).to_edge(DOWN)

        self.play(
            FadeOut(VGroup(formula, split_note, bar1, bar2, b1_label, b2_label, growth), shift=UP * 0.2),
            Transform(title, pivot_title),
            run_time=0.8,
        )
        self.wait(3)
        self.play(FadeIn(q, shift=UP * 0.12), run_time=0.7)
        self.wait(3)
        self.play(FadeIn(t1, shift=UP * 0.1), Write(t1_t), run_time=0.6)
        self.play(FadeIn(t2, shift=UP * 0.1), Write(t2_t), run_time=0.6)
        self.play(FadeIn(t3, shift=UP * 0.1), Write(t3_t), run_time=0.6)
        self.play(FadeIn(close, shift=UP * 0.1), run_time=0.7)
        self.wait(4)

class Chapter2Solution(Scene):
    def construct(self):
        title_page = Text("Chapter 2: Solution", font_size=48)
        subtitle = Paragraph(
            "Memory-Bounded Exact Search\nPython Control Plane + C++ Compute Plane",
            font_size=26,
            alignment="center",
        ).next_to(title_page, DOWN, buff=0.6)
        cover = VGroup(title_page, subtitle).center()
        self.play(Write(cover), run_time=1.2)
        self.wait(2.5)
        self.play(FadeOut(cover), run_time=0.6)

        title = Text("k = 8 Changes the Problem", font_size=42).to_edge(UP)
        self.play(Write(title), run_time=0.9)

        left = RoundedRectangle(width=5.6, height=2.5, corner_radius=0.15, color=RED_C).shift(LEFT * 3.4 + DOWN * 0.4)
        left_t = Text("Classic 4+4 MITM", font_size=28, color=RED_A).move_to(left.get_top() + DOWN * 0.3)
        left_b = Paragraph(
            "Precompute all 4-term sums\nTable size ~ O(MaxVal^4)",
            font_size=24,
            alignment="center",
            line_spacing=0.8,
        ).move_to(left.get_center() + DOWN * 0.1)

        right = RoundedRectangle(width=5.6, height=2.5, corner_radius=0.15, color=YELLOW_D).shift(RIGHT * 3.4 + DOWN * 0.4)
        right_t = Text("Observed Bottleneck", font_size=28, color=YELLOW_A).move_to(right.get_top() + DOWN * 0.3)
        right_b = Paragraph(
            "RAM explodes before\nmath becomes interesting",
            font_size=24,
            alignment="center",
            line_spacing=0.8,
        ).move_to(right.get_center() + DOWN * 0.1)

        link = Arrow(left.get_right(), right.get_left(), buff=0.2, color=GRAY_B, stroke_width=5)
        self.play(FadeIn(VGroup(left, left_t, left_b), shift=RIGHT * 0.2), run_time=0.8)
        self.play(FadeIn(VGroup(right, right_t, right_b), shift=LEFT * 0.2), GrowArrow(link), run_time=0.8)
        self.wait(3)

        py_cpp_title = Text("Optimization Strategy", font_size=40, color=YELLOW_B).to_edge(UP)
        py_box = RoundedRectangle(width=5.2, height=2.9, corner_radius=0.15, color=BLUE_C).shift(LEFT * 3.2 + DOWN * 0.3)
        py_head = Text("Python Control Plane", font_size=28, color=BLUE_A).move_to(py_box.get_top() + DOWN * 0.3)
        py_body = Paragraph(
            "UI + parameters\nexperiment logging\nharness orchestration",
            font_size=23,
            alignment="center",
            line_spacing=0.8,
        ).move_to(py_box.get_center() + DOWN * 0.15)

        cpp_box = RoundedRectangle(width=5.2, height=2.9, corner_radius=0.15, color=TEAL_C).shift(RIGHT * 3.2 + DOWN * 0.3)
        cpp_head = Text("C++ Compute Plane", font_size=28, color=TEAL_A).move_to(cpp_box.get_top() + DOWN * 0.3)
        cpp_body = Paragraph(
            "Hot loop for m=8\nsearch_m8_blocked_sieved\nblocked scan + modular filters",
            font_size=23,
            alignment="center",
            line_spacing=0.8,
        ).move_to(cpp_box.get_center() + DOWN * 0.15)

        callout = Text("Exact solver stays exact; Python overhead leaves critical path", font_size=24, color=GRAY_A).to_edge(DOWN)
        bridge = Arrow(py_box.get_right(), cpp_box.get_left(), buff=0.2, color=GRAY_B, stroke_width=5)

        self.play(FadeOut(VGroup(left, left_t, left_b, right, right_t, right_b, link), shift=UP * 0.2), Transform(title, py_cpp_title), run_time=0.8)
        self.play(FadeIn(VGroup(py_box, py_head, py_body), shift=RIGHT * 0.2), run_time=0.7)
        self.play(FadeIn(VGroup(cpp_box, cpp_head, cpp_body), shift=LEFT * 0.2), GrowArrow(bridge), run_time=0.7)
        self.play(FadeIn(callout, shift=UP * 0.1), run_time=0.6)
        self.wait(3)

        algo_title = Text("Index + Cheap Filters + Exact Check", font_size=38, color=YELLOW_B).to_edge(UP)
        step1 = RoundedRectangle(width=3.9, height=1, corner_radius=0.12, color=GREEN_C).shift(UP * 1.85)
        step1_t = Text("Build 4-sum index", font_size=22).move_to(step1)
        step2 = RoundedRectangle(width=3.9, height=1, corner_radius=0.12, color=GREEN_C).shift(UP * 0.6)
        step2_t = Text("Bucket by mod 2^16", font_size=22).move_to(step2)
        step3 = RoundedRectangle(width=3.9, height=1, corner_radius=0.12, color=GREEN_C).shift(DOWN * 0.65)
        step3_t = Text("Attach two prime residues", font_size=22).move_to(step3)
        step4 = RoundedRectangle(width=3.9, height=1, corner_radius=0.12, color=GREEN_C).shift(DOWN * 1.9)
        step4_t = Text("Exact equality verify", font_size=22).move_to(step4)
        f1 = Arrow(step1.get_bottom(), step2.get_top(), buff=0.15, color=GRAY_B)
        f2 = Arrow(step2.get_bottom(), step3.get_top(), buff=0.15, color=GRAY_B)
        f3 = Arrow(step3.get_bottom(), step4.get_top(), buff=0.15, color=GRAY_B)
        receipts = Text("Like matching two receipts to hit one total", font_size=25, color=GRAY_A).to_edge(DOWN)

        self.play(FadeOut(VGroup(py_box, py_head, py_body, cpp_box, cpp_head, cpp_body, bridge, callout), shift=UP * 0.2), Transform(title, algo_title), run_time=0.8)
        self.play(FadeIn(step1), Write(step1_t), run_time=0.6)
        self.play(GrowArrow(f1), FadeIn(step2), Write(step2_t), run_time=0.6)
        self.play(GrowArrow(f2), FadeIn(step3), Write(step3_t), run_time=0.6)
        self.play(GrowArrow(f3), FadeIn(step4), Write(step4_t), run_time=0.6)
        self.play(FadeIn(receipts, shift=UP * 0.12), run_time=0.6)
        self.wait(3)

        val_title = Text("Validation Culture", font_size=40, color=YELLOW_B).to_edge(UP)
        v1 = RoundedRectangle(width=8.8, height=1.15, corner_radius=0.12, color=BLUE_C).shift(UP * 1.4)
        v1_t = Text("Small bounds: brute vs Python MITM must match", font_size=25).move_to(v1)
        v2 = RoundedRectangle(width=8.8, height=1.15, corner_radius=0.12, color=BLUE_C).shift(UP * 0.1)
        v2_t = Text("C++ enabled: brute vs Python MITM vs C++ must match", font_size=25).move_to(v2)
        v3 = RoundedRectangle(width=8.8, height=1.15, corner_radius=0.12, color=BLUE_C).shift(DOWN * 1.2)
        v3_t = Text("Fail loud on mismatch; single clear success signal", font_size=25).move_to(v3)

        self.play(FadeOut(VGroup(step1, step1_t, step2, step2_t, step3, step3_t, step4, step4_t, f1, f2, f3, receipts), shift=UP * 0.2), Transform(title, val_title), run_time=0.8)
        self.play(FadeIn(v1, shift=UP * 0.1), Write(v1_t), run_time=0.6)
        self.play(FadeIn(v2, shift=UP * 0.1), Write(v2_t), run_time=0.6)
        self.play(FadeIn(v3, shift=UP * 0.1), Write(v3_t), run_time=0.6)
        self.wait(3)

        fix_title = Text("Prototype Fixes -> Final Correctness", font_size=36, color=YELLOW_B).to_edge(UP)
        fix1 = Text("Old prototypes: close structure, edge-case bugs", font_size=24, color=WHITE).shift(UP * 1.4)
        fix2 = Text("Primitive grouping fix: divide BOTH left terms and b by global gcd", font_size=24, color=GREEN_C).shift(UP * 0.5)
        fix3 = Text("Interface cleanup: callback/progress mismatches resolved", font_size=24, color=WHITE).shift(DOWN * 0.4)
        fix4 = Text("Correctness no longer depends on post-processing", font_size=24, color=GREEN_C).shift(DOWN * 1.3)

        self.play(FadeOut(VGroup(v1, v1_t, v2, v2_t, v3, v3_t), shift=UP * 0.2), Transform(title, fix_title), run_time=0.8)
        self.play(Write(fix1), run_time=0.6)
        self.play(Write(fix2), run_time=0.6)
        self.play(Write(fix3), run_time=0.6)
        self.play(Write(fix4), run_time=0.6)
        self.wait(3)

        out_title = Text("Chapter 2 Outcome", font_size=40, color=YELLOW_B).to_edge(UP)
        p1 = RoundedRectangle(width=4.0, height=1.25, corner_radius=0.12, color=TEAL_C).shift(UP * 2)
        p1_t = Text("Index vs Query", font_size=24).move_to(p1)
        p2 = RoundedRectangle(width=4.0, height=1.25, corner_radius=0.12, color=TEAL_C)
        p2_t = Text("Cheap before Expensive", font_size=24).move_to(p2)
        p3 = RoundedRectangle(width=4.0, height=1.25, corner_radius=0.12, color=TEAL_C).shift(DOWN * 2)
        p3_t = Text("Fast + Correct + Scalable", font_size=24).move_to(p3)
        final_line = Text(
            "Same systems instinct used in retrieval, training, and inference infrastructure.",
            font_size=24,
            color=GRAY_A,
        ).to_edge(DOWN)

        self.play(FadeOut(VGroup(fix1, fix2, fix3, fix4), shift=UP * 0.2), Transform(title, out_title), run_time=0.8)
        self.play(FadeIn(p1, shift=UP * 0.1), Write(p1_t), run_time=0.6)
        self.play(FadeIn(p2, shift=UP * 0.1), Write(p2_t), run_time=0.6)
        self.play(FadeIn(p3, shift=UP * 0.1), Write(p3_t), run_time=0.6)
        self.play(FadeIn(final_line, shift=UP * 0.1), run_time=0.6)
        self.wait(4)

class Chapter3Applications(Scene):
    def construct(self):
        def stage_box(label, width=3.6, height=1.1, color=BLUE_C, font_size=22):
            box = RoundedRectangle(width=width, height=height, corner_radius=0.12, color=color)
            txt = Tex(label, font_size=font_size, tex_environment="center").move_to(box)
            return VGroup(box, txt)

        chapter_title = Text("Chapter 3: Applications", font_size=46)
        subtitle = Paragraph(
            "Math Research and the Use in AI Systems Development",
            font_size=26,
            alignment="center",
        ).next_to(chapter_title, DOWN, buff=0.6)
        title_page = VGroup(chapter_title, subtitle)

        self.play(Write(title_page), run_time=1.0)
        self.wait(5)

        # Paragraph 1: AI tradeoffs + exactness contract
        p1_title = Text("Same Tradeoffs as AI Systems", font_size=36, color=YELLOW_B).to_edge(UP)
        p1_left = stage_box("Precompute Index", color=TEAL_C, font_size=30).shift(LEFT * 3.8)
        p1_right = stage_box("Stream Compute", color=BLUE_C, font_size=30).shift(RIGHT * 3.8)
        p1_arrow = Arrow(p1_left.get_right(), p1_right.get_left(), buff=0.15, color=GRAY_B)
        p1_m1 = stage_box("Throughput vs Memory", width=4.4, color=GREEN_C, font_size=30).next_to(p1_title, DOWN)
        p1_badge = RoundedRectangle(width=12.5, height=1.2, corner_radius=0.14, color=YELLOW_D)
        p1_badge_t = Text(
            "Exact integer arithmetic only: no floating-point heuristics, no ad hoc constraints",
            font_size=22,
            color=YELLOW_A,
        ).move_to(p1_badge)
        p1_badge_whole = VGroup(p1_badge_t, p1_badge).to_edge(DOWN)
        p1_m2 = stage_box("Optimize Without Losing Correctness", width=6.1, color=GREEN_C, font_size=30).to_edge(DOWN)

        self.play(ReplacementTransform(title_page, p1_title), run_time=0.6)
        title_page = p1_title
        self.play(FadeIn(p1_left, shift=RIGHT * 0.2), FadeIn(p1_right, shift=LEFT * 0.2), GrowArrow(p1_arrow), run_time=0.8)
        self.play(FadeIn(p1_m1, shift=UP * 0.12), FadeIn(p1_m2, shift=UP * 0.12), run_time=0.7)
        self.wait(5)
        self.play(Transform(p1_m2, p1_badge_whole), run_time=0.7)
        self.wait(3)

        # Paragraph 2: MITM as database indexing
        p2_title = Text("MITM as Database Indexing", font_size=36, color=YELLOW_B).to_edge(UP)
        left_half = stage_box("Half A: precompute \n nondecreasing r-tuples", width=4, color=TEAL_C, font_size=20).shift(LEFT * 4 + DOWN * 0.2)
        index_tbl = RoundedRectangle(width=2, height=2, corner_radius=0.16, color=YELLOW_D).shift(DOWN * 0.2)
        index_txt = Tex("Index\nTable", font_size=28, color=YELLOW_A, tex_environment="center").move_to(index_tbl)
        right_half = stage_box("Half B: query \n complementary tuples \n per b", width=4, color=BLUE_C, font_size=20).shift(RIGHT * 4 + DOWN * 0.2)
        a1 = Arrow(left_half.get_right(), index_tbl.get_left(), buff=0.15, color=GRAY_B)
        a2 = Arrow(right_half.get_left(), index_tbl.get_right(), buff=0.15, color=GRAY_B)
        p2_note = Text(
            "Build once, query many times: same structure as retrieval indexing",
            font_size=24,
            color=GRAY_A,
        ).to_edge(DOWN)

        self.play(ReplacementTransform(title_page, p2_title))
        title_page = p2_title
        self.play(
            FadeOut(VGroup(p1_left, p1_right, p1_arrow, p1_m1, p1_m2, p1_badge, p1_badge_t), shift=UP * 0.2),
            run_time=0.8,
        )
        self.play(FadeIn(left_half, shift=RIGHT * 0.2), FadeIn(index_tbl), Write(index_txt), run_time=0.8)
        self.play(GrowArrow(a1), FadeIn(right_half, shift=LEFT * 0.2), GrowArrow(a2), run_time=0.8)
        self.play(FadeIn(p2_note, shift=UP * 0.1), run_time=0.6)
        self.wait(3)

        # Paragraph 3: cheap filter then expensive check (RAG analogy + bitmask sieve)
        p3_title = Text("Cheap Filter -> Expensive Check", font_size=36, color=YELLOW_B).to_edge(UP)
        rag1 = stage_box("RAG: fast retrieval", width=3.0, color=PURPLE_C, font_size=25).shift(LEFT * 4.4 + UP * 1.3)
        rag2 = stage_box("RAG: rerank", width=3.0, color=PURPLE_C, font_size=25).shift(UP * 1.3)
        rag3 = stage_box("RAG: exact score", width=3.0, color=PURPLE_C, font_size=25).shift(RIGHT * 4.4 + UP * 1.3)
        r12 = Arrow(rag1.get_right(), rag2.get_left(), buff=0.12, color=GRAY_B)
        r23 = Arrow(rag2.get_right(), rag3.get_left(), buff=0.12, color=GRAY_B)

        sv1 = stage_box("Solver: residue sieve", width=3.0, color=GREEN_C, font_size=25).shift(LEFT * 4.4 + DOWN * 0.5)
        sv2 = stage_box("Solver: candidate set", width=3.0, color=GREEN_C, font_size=25).shift(DOWN * 0.5)
        sv3 = stage_box("Solver: exact verify", width=3.0, color=GREEN_C, font_size=25).shift(RIGHT * 4.4 + DOWN * 0.5)
        s12 = Arrow(sv1.get_right(), sv2.get_left(), buff=0.12, color=GRAY_B)
        s23 = Arrow(sv2.get_right(), sv3.get_left(), buff=0.12, color=GRAY_B)

        bits = VGroup(*[Square(side_length=0.28, color=BLUE_C) for _ in range(12)]).arrange(RIGHT, buff=0.04).to_edge(DOWN).shift(UP * 0.8)
        on_idx = {1, 3, 4, 8, 10}
        for i, sq in enumerate(bits):
            if i in on_idx:
                sq.set_fill(BLUE_C, opacity=0.85)
            else:
                sq.set_fill(BLACK, opacity=0.0)
        bits_lbl = Text("Bitmask residue sumset (per prime p)", font_size=22, color=GRAY_A).next_to(bits, UP, buff=0.18)
        impossible = MathTex(r"\text{Reject impossible } b^k mod p \text{ early}", font_size=22, color=RED_C).next_to(bits, DOWN, buff=0.2)

        self.play(ReplacementTransform(title_page, p3_title))
        title_page = p3_title
        self.play(
            FadeOut(VGroup(left_half, index_tbl, index_txt, right_half, a1, a2, p2_note), shift=UP * 0.2),
            run_time=0.8,
        )
        self.play(FadeIn(rag1), GrowArrow(r12), FadeIn(rag2), GrowArrow(r23), FadeIn(rag3), run_time=0.8)
        self.play(FadeIn(sv1), GrowArrow(s12), FadeIn(sv2), GrowArrow(s23), FadeIn(sv3), run_time=0.8)
        self.play(FadeIn(bits_lbl, shift=UP * 0.1), FadeIn(bits), FadeIn(impossible, shift=UP * 0.1), run_time=0.8)
        self.wait(3)

        # Paragraph 4: correctness harness discipline
        p4_title = Text("Performance Only Counts If Correct", font_size=36, color=YELLOW_B).to_edge(UP)
        ref = stage_box("Reference: brute force", width=3.25, color=BLUE_C, font_size=25).shift(LEFT * 4 + UP * 0.7)
        py = stage_box("Python MITM", width=3.25, color=TEAL_C, font_size=25).shift(UP * 0.7)
        cpp = stage_box("C++ optimized path", width=3.25, color=GREEN_C, font_size=25).shift(RIGHT * 4 + UP * 0.7)
        c1 = Arrow(ref.get_right(), py.get_left(), buff=0.12, color=GRAY_B)
        c2 = Arrow(py.get_right(), cpp.get_left(), buff=0.12, color=GRAY_B)
        gate = RoundedRectangle(width=8.6, height=1.4, corner_radius=0.14, color=YELLOW_D).shift(DOWN * 1.0)
        gate_t = Text("Correctness gate: all paths must match known cases", font_size=24, color=YELLOW_A).move_to(gate)
        fail_loud = Text("Any mismatch fails loudly before scaling", font_size=24, color=RED_C).to_edge(DOWN)

        self.play(ReplacementTransform(title_page, p4_title))
        title_page = p4_title
        self.play(
            FadeOut(VGroup(rag1, rag2, rag3, r12, r23, sv1, sv2, sv3, s12, s23, bits, bits_lbl, impossible), shift=UP * 0.2),
            run_time=0.8,
        )
        self.play(FadeIn(ref), FadeIn(py), FadeIn(cpp), GrowArrow(c1), GrowArrow(c2), run_time=0.8)
        self.play(FadeIn(gate, shift=UP * 0.1), Write(gate_t), run_time=0.7)
        self.play(FadeIn(fail_loud, shift=UP * 0.1), run_time=0.6)
        self.wait(3)

        # Paragraph 5: memory bandwidth bottleneck + transferable techniques
        p5_title = Text("Scaling Barrier: Memory Bandwidth", font_size=36, color=YELLOW_B).to_edge(UP)
        comp_bar = Rectangle(width=1.8, height=1.4, color=BLUE_C, fill_opacity=0.7).shift(LEFT * 3.2 + DOWN * 2.2)
        mem_bar = Rectangle(width=1.8, height=3.8, color=RED_C, fill_opacity=0.7).shift(LEFT * 0.9 + UP * -1.0)
        comp_lbl = Text("Compute", font_size=24).next_to(comp_bar, DOWN, buff=0.2)
        mem_lbl = Text("Memory move", font_size=24).next_to(mem_bar, DOWN, buff=0.2)
        ai_map = Paragraph(
            "AI analogs:\nactivations, KV cache,\noptimizer states",
            font_size=23,
            alignment="left",
            line_spacing=0.85,
            color=GRAY_A,
        ).shift(RIGHT * 3.6 + UP * 2)
        tech1 = stage_box("Bucketing", width=3.2, color=TEAL_C, font_size=30).shift(RIGHT * 3.6 + DOWN * 0.1)
        tech2 = stage_box("Compact reps", width=3.2, color=TEAL_C, font_size=30).shift(RIGHT * 3.6 + DOWN * 1.3)
        tech3 = stage_box("Blocked chunks", width=3.2, color=TEAL_C, font_size=30).shift(RIGHT * 3.6 + DOWN * 2.5)

        self.play(ReplacementTransform(title_page, p5_title))
        title_page = p5_title
        self.play(
            FadeOut(VGroup(ref, py, cpp, c1, c2, gate, gate_t, fail_loud, p4_title), shift=UP * 0.2),
            run_time=0.8,
        )
        self.play(FadeIn(comp_bar), FadeIn(mem_bar), FadeIn(comp_lbl), FadeIn(mem_lbl), run_time=0.7)
        self.play(FadeIn(ai_map, shift=UP * 0.1), run_time=0.6)
        self.play(FadeIn(tech1, shift=UP * 0.1), FadeIn(tech2, shift=UP * 0.1), FadeIn(tech3, shift=UP * 0.1), run_time=0.8)
        self.wait(3)

        # Paragraph 6: startup workflow and measurable metrics
        p6_title = Text("Applied Startup Workflow", font_size=36, color=YELLOW_B).to_edge(UP)
        w1 = stage_box("1) Find bottleneck", width=3.2, color=BLUE_C, font_size=25).shift(LEFT * 3 + UP * 1.3)
        w2 = stage_box("2) Make measurable", width=3.2, color=BLUE_C, font_size=25).shift(RIGHT * 3 + UP * 1.3)
        w3 = stage_box("3) Redesign for constraints", width=3.2, color=BLUE_C, font_size=25).shift(RIGHT * 3 + DOWN * 0.8)
        w4 = stage_box("4) Prove correctness", width=3.2, color=BLUE_C, font_size=25).shift(LEFT * 3 + DOWN * 0.8)
        wa1 = Arrow(w1.get_right(), w2.get_left(), buff=0.12, color=GRAY_B)
        wa2 = Arrow(w2.get_bottom(), w3.get_top(), buff=0.12, color=GRAY_B)
        wa3 = Arrow(w3.get_left(), w4.get_right(), buff=0.12, color=GRAY_B)
        targets = Text("Instrument: embedding index, classifier pipeline, scheduler, or document workflow", font_size=22, color=GRAY_A).shift(DOWN * 1.8)
        metrics = RoundedRectangle(width=9.5, height=1.25, corner_radius=0.14, color=GREEN_C).to_edge(DOWN)
        metrics_t = Text("User-critical metrics: latency and cost (memory/compute)", font_size=24).move_to(metrics)

        self.play(ReplacementTransform(title_page, p6_title))
        title_page = p6_title
        self.play(
            FadeOut(VGroup(comp_bar, mem_bar, comp_lbl, mem_lbl, ai_map, tech1, tech2, tech3), shift=UP * 0.2),
            run_time=0.8,
        )
        self.play(FadeIn(w1), GrowArrow(wa1), FadeIn(w2), GrowArrow(wa2), FadeIn(w3), GrowArrow(wa3), FadeIn(w4), run_time=0.9)
        self.wait(5)
        self.play(FadeIn(targets, shift=UP * 0.1), run_time=0.6)
        self.play(FadeIn(metrics, shift=UP * 0.1), Write(metrics_t), run_time=0.6)
        self.wait(3)

        # Paragraph 7: week-one deliverable memo + patch + test
        p7_title = Text("Week One Delivery", font_size=36, color=YELLOW_B).to_edge(UP)
        memo = RoundedRectangle(width=4.75, height=3, corner_radius=0.14, color=YELLOW_D).shift(UP)
        memo_title = Text("Engineering Memo + Patch", font_size=25, color=YELLOW_A).move_to(memo.get_top() + DOWN * 0.25)
        lines = Paragraph("1) Current behavior "
                            "\n 2) Bottleneck " \
                            "\n 3) Hypothesis " \
                            "\n 4) Smallest Change " \
                            "\n 5) Test: Optimized Path",
                            line_spacing=0.7, 
                            alignment="center",
                            font_size=20).move_to(memo.get_center()).shift(DOWN * 0.2) 
        
        ship = stage_box("It only ships if it passes the correctness gate", width=5.5, color=GREEN_C, font_size=25).shift(DOWN * 2)
        ship_list = Paragraph(
            "Examples: new kernels caching layers rerankers quantization paths",
            font_size=21,
            alignment="left",
            line_spacing=0.82,
            color=GRAY_A,
        ).to_edge(DOWN)

        self.play(ReplacementTransform(title_page, p7_title))
        title_page = p7_title
        self.play(
            FadeOut(VGroup(w1, w2, w3, w4, wa1, wa2, wa3, targets, metrics, metrics_t), shift=UP * 0.2),
            run_time=0.8,
        )
        self.play(ReplacementTransform(title_page, p7_title))
        self.play(FadeIn(memo, shift=UP * 0.1), Write(memo_title), run_time=0.7)
        self.play(Write(lines), run_time=1.0)
        self.play(FadeIn(ship, shift=UP * 0.1), FadeIn(ship_list, shift=UP * 0.1), run_time=0.7)
        self.wait(3)

        # Paragraph 8: transfer patterns across AI systems
        p8_title = Text("Transferable Candidate-Reduction Pattern", font_size=36, color=YELLOW_B).to_edge(UP)
        col1_a = stage_box("Embeddings: ANN retrieval", width=3.4, color=TEAL_C, font_size=20).shift(LEFT * 4.5 + UP * 1.75)
        col1_b = stage_box("Rerank", width=3.2, color=TEAL_C, font_size=20).next_to(col1_a, DOWN, buff=0.5)
        col1_c = stage_box("Exact scoring", width=3, color=TEAL_C, font_size=20).next_to(col1_b, DOWN, buff=0.5)
        c1a = Arrow(col1_a.get_bottom(), col1_b.get_top(), buff=0.07, color=GRAY_B)
        c1b = Arrow(col1_b.get_bottom(), col1_c.get_top(), buff=0.07, color=GRAY_B)

        col2_a = stage_box("Workflow: rule filter", width=3.4, color=BLUE_C, font_size=20).shift(UP * 1.75)
        col2_b = stage_box("Model classify", width=3.2, color=BLUE_C, font_size=20).next_to(col2_a, DOWN, buff=0.5)
        col2_c = stage_box("Human review edge cases", width=3, color=BLUE_C, font_size=20).next_to(col2_b, DOWN, buff=0.5)
        c2a = Arrow(col2_a.get_bottom(), col2_b.get_top(), buff=0.07, color=GRAY_B)
        c2b = Arrow(col2_b.get_bottom(), col2_c.get_top(), buff=0.07, color=GRAY_B)

        col3_a = stage_box("Scheduling: prune constraints", width=3.4, color=PURPLE_C, font_size=20).shift(RIGHT * 4.5 + UP * 1.75)
        col3_b = stage_box("Generate candidates", width=3.2, color=PURPLE_C, font_size=20).next_to(col3_a, DOWN, buff=0.5)
        col3_c = stage_box("Final verification", width=3.0, color=PURPLE_C, font_size=20).next_to(col3_b, DOWN, buff=0.5)
        c3a = Arrow(col3_a.get_bottom(), col3_b.get_top(), buff=0.07, color=GRAY_B)
        c3b = Arrow(col3_b.get_bottom(), col3_c.get_top(), buff=0.07, color=GRAY_B)
        p8_footer = Text("Aggressively shrink candidates, never break correctness contract", font_size=24, color=GRAY_A).to_edge(DOWN)

        self.play(ReplacementTransform(title_page, p8_title))
        title_page = p8_title
        self.play(
            FadeOut(VGroup(memo, memo_title, lines, ship, ship_list), shift=UP * 0.2),
            run_time=0.8,
        )
        self.play(FadeIn(col1_a), FadeIn(col2_a), FadeIn(col3_a), run_time=0.6)
        self.play(GrowArrow(c1a), FadeIn(col1_b), GrowArrow(c2a), FadeIn(col2_b), GrowArrow(c3a), FadeIn(col3_b), run_time=0.7)
        self.play(GrowArrow(c1b), FadeIn(col1_c), GrowArrow(c2b), FadeIn(col2_c), GrowArrow(c3b), FadeIn(col3_c), run_time=0.7)
        self.play(FadeIn(p8_footer, shift=UP * 0.1), run_time=0.6)
        self.wait(3)

        # Paragraph 9: system hygiene and trustworthy production behavior
        p9_title = Text("System Hygiene for Production AI", font_size=36, color=YELLOW_B).to_edge(UP)
        h1 = stage_box("Progress reporting", width=3.25, color=GREEN_C).shift(LEFT * 5 + UP * 0.8)
        h2 = stage_box("Reproducible runs", width=3.25, color=GREEN_C).shift(UP * 0.8)
        h3 = stage_box("Clear failure modes", width=3.25, color=GREEN_C).shift(RIGHT * 5 + UP * 0.8)
        h4 = stage_box("Observable and testable", width=4, color=GREEN_C).shift(LEFT * 2.5 + DOWN * 1.5)
        h5 = stage_box("Safe to optimize and change", width=4, color=GREEN_C).shift(RIGHT * 2.5 + DOWN * 1.5)
        final_msg = Text(
            "Make it fast, make it transparent, make it trustworthy.",
            font_size=28,
            color=YELLOW_A,
        ).to_edge(DOWN)

        self.play(ReplacementTransform(title_page, p9_title))
        title_page = p9_title
        self.play(
            FadeOut(VGroup(col1_a, col1_b, col1_c, col2_a, col2_b, col2_c, col3_a, col3_b, col3_c, c1a, c1b, c2a, c2b, c3a, c3b, p8_footer), shift=UP * 0.2),
            run_time=0.8,
        )
        self.play(FadeIn(h1, shift=UP * 0.1), FadeIn(h2, shift=UP * 0.1), FadeIn(h3, shift=UP * 0.1), run_time=0.7)
        self.play(FadeIn(h4, shift=UP * 0.1), FadeIn(h5, shift=UP * 0.1), run_time=0.7)
        self.play(FadeIn(final_msg, shift=UP * 0.1), run_time=0.7)
        self.wait(4)