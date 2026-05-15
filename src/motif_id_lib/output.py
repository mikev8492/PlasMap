
import sys, re, math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# --- Module-level constants ------------------
import random as _random


# Create random color palette for enzymes:
# 1. create fixed random seed: 
_random.seed(42)
# 2. take all colours from tab20, tab20b, tab20c (60 total)
_tab20 = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors)
# 3. shuffle them with a fixed random seed
_random.shuffle(_tab20)
# assign one distinct colour per enzyme
PALETTE = _tab20

# Ambiguity codes:
IUPAC = {
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}

# Set hex colors to each base:
BASE_COLORS  = {'A': '#69db7c', 'T': '#ff6b6b', 'G': '#ffd43b', 'C': '#4dabf7'}

# Complement base pairing:
COMPLEMENT   = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

# sys.stdout colors:
HIGHLIGHT = "\033[36m" # teal
RESET = "\033[0m"

# ===============================
# Shared data-preparation mixin
# ===============================
class _EnzymeDataMixin:
    """
    Purpose:
    --------
        Helper class for building:
            - colour maps 
            - motif index lookup table
            - cut-position (top and bottom) tables 
        
        Inherits the following from parent class:
            self.results, self.seq, self.seq_len, self.enzyme_names
    """

    def _build_color_map(self) -> dict:
        """
        Purpose:
        --------
            Assign a palette colour to each enzyme name.

        Returns:
        --------
            {enzyme name: (r, g, b)}

        """
        return {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(self.enzyme_names)}

    # NOTE:@staticmethod allows function to be standalone (doesnt need "self" argument)
    @staticmethod  
    def _top_cut_offset(cut_notation: str) -> int:
        """
        Purpose:
        --------
            Offset of the top-strand cut within the motif (^ position, _ stripped).

        Arguments:
        ----------
            cut_notation : str
                Enzyme cut notation string containing '^' for the top-strand cut position and '_' for the bottom-strand cut position.
                - Example: 'G^AATTC_'

        Returns:
        --------
            offset index: integer
        """
        return cut_notation.replace('_', '').index('^')

    @staticmethod
    def _bot_cut_offset(cut_notation: str) -> int:
        """
        Purpose:
        --------
            Offset of the bottom-strand cut within the motif (_ position, ^ stripped).

        Arguments:
        ----------
            cut_notation : str

        Returns:
        --------
            offset index: integer
        """
        return cut_notation.replace('^', '').index('_')

    def _build_char_enzyme_map(self) -> dict:
        """
        Purpose:
        --------
            Builds lookup table with every sequence index that falls within the enzyme motif using:
                Results start index + motif length
        Returns:
        --------
            {index: enzyme name}
        """
        char_enzyme = {}
        for name, data in self.results.items():
            motif_len = len(data[0])
            for site_pos in data[3]:
                for i in range(site_pos, site_pos + motif_len):
                    char_enzyme[i] = name
        return char_enzyme

    def _build_top_cut_positions(self) -> dict:
        """
        Purpose:
        --------
            Build lookup table for top strand cuts:
                absolute cut position = site + top offset

        Returns:
        --------
            Lookup table:
                {top cut index: [enzyme names]}
                Dictionary mapping each absolute top-strand cut position to a list
                of enzyme names that cut at that position. Enzymes with zero cuts
                are excluded
        """
        positions: dict[int, list[str]] = {}
        for name, data in self.results.items():
            _, notation, cut_count, sites = data # unpacks data iterable(tuple)
            if cut_count == 0:
                continue
            offset = self._top_cut_offset(notation)
            for site in sites:
                positions.setdefault(site + offset, []).append(name)
        return positions

    def _build_bot_cut_positions(self) -> dict:
        """
        Purpose:
        --------
            Build lookup table for top strand cuts:
                absolute cut position = site + bottom offset
        Returns:
        --------
            Lookup table:
                {bottom cut index: [enzyme names]}
        """
        positions: dict[int, list[str]] = {}
        for name, data in self.results.items():
            _, notation, cut_count, sites = data
            if cut_count == 0:
                continue
            offset = self._bot_cut_offset(notation)
            for site in sites:
                positions.setdefault(site + offset, []).append(name)
        return positions

    def _build_legend_handles(self, color_map: dict) -> list:
        """
        Purpose:
        --------
            Builds handles: 
                label and Patch (colored box) for enzyme legend, used by matplotlib
        
        Arguments:
        ----------
            color_map : dict
                {enzyme name: (r, g, b)} mapping returned by _build_color_map().
        
        Returns:
        --------
            List of matplotlib Patch objects, one per cutting enzyme.
            Each Patch carries:
                - facecolor / edgecolor matching the enzyme's assigned palette colour
                - label formatted as: "EcoRI [GAATTC] @ 283, 284, 300"
            Enzymes with zero cuts are excluded.
        """
        handles = []
        for name, data in self.results.items():
            if data[2] == 0:
                continue
            pos_str = ', '.join(str(p) for p in data[3])
            label   = f"{name} [{data[0]}] @ {pos_str}"
            handles.append(mpatches.Patch(
                facecolor=color_map[name], edgecolor=color_map[name],
                alpha=0.8, label=label,
            ))
        return handles

# ==================
# CircularMap
# ==================
class CircularMap(_EnzymeDataMixin):
    """
    Purpose:
    --------
        Generate Circular plasmid map with restriction-enzyme cut sites.

    Parameters
    ----------
    results : dict
        { enzyme_name: [recognition_seq, cut_notation, cut_count, [positions]] }
    plasmid_sequence : str
        Full nucleotide sequence.
    title : str
        Text printed at the centre of the map.
    figsize : tuple
        Figure dimensions in inches.
    dpi : int
        Output resolution.
    output_path : str
        Destination file for the saved figure.
    """

    # --- Geometry constants ------------------------------
    RADIUS    : float = 1.6    # backbone circle radius
    TICK_R    : float = 0.12   # radial length of cut-site tick marks
    LBL_GAP   : float = 0.22   # extra radial gap from tick to label
    MC_GAP    : float = 0.07   # extra radial gap for multi-cut enzymes
    BB_LW     : int   = 3      # backbone linewidth
    POS_TICK_INTERVAL: int = 500  # bp between ruler ticks on the backbone

    def __init__(self, results: dict, plasmid_sequence: str, title: str = "Plasmid", figsize: tuple = (10, 10), dpi: int = 150, output_path: str = "results/circular_map.png") -> None:
        self.results      = results
        self.seq          = plasmid_sequence.upper().replace(' ', '').replace('\n', '')
        self.seq_len      = len(self.seq)
        self.title        = title
        self.figsize      = figsize
        self.DPI          = dpi
        self.OUTPUT_PATH  = output_path

        self.enzyme_names = list(results.keys())
        self.color_map    = self._build_color_map()
        self.cutters      = {k: v for k, v in results.items() if v[2] > 0} # only non-zero cuts

    # --- Figure setup -------------------------------------------------

    def _make_figure(self) -> tuple:
        """
        Purpose:
        --------
            Create base figure for the plot:
                - Create square canvas (self.figsize defaults to 10x10)
                - set aspect to equal (same ratio for x and y) to prevent stretching the circle
                - Turn off default matplotlib axes features (border box, tick marks and labels, etc.)
                - Set axes limits outside backbone circle to allow enzyme labels
                - Set figure background
        Returns:
        --------
            Two objects to build one panel in the plot:
            - (fig: figure object, ax: axes object)
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_xlim(-3.2, 3.2)
        ax.set_ylim(-3.0, 3.0)
        fig.patch.set_facecolor('white')
        return fig, ax

    # --- Drawing helpers ---------------------------------------------

    def _draw_backbone(self, ax) -> None:
        """
        Purpose:
        --------
            Add an open circle patch to the axes with set radius and line width constants. 
            (patches are gathered and drawn at render)

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The axes object on which the backbone circle is drawn.
        """
        ax.add_patch(plt.Circle(
            (0, 0), self.RADIUS,
            fill=False, color='#444', linewidth=self.BB_LW,
        ))

    def _draw_position_ticks(self, ax) -> None:
        """
        Purpose:
        --------
            Draw evenly spaced bp-position ticks and labels around the backbone.
            - Convert base pair positions to angles on the circle (convert to radians for np.cos/np.sin)
                - bp position 0 = 0 degrees, 1/2 sequence length = 180 degrees, last bp position = 360 degrees, etc.
            - draw tick mark at calculated angles (from inner radius r0 to outer radius r1)
            - draw bp number outside tick mark (r1 + 0.7)
        """
        R = self.RADIUS
        for bp in range(0, self.seq_len, self.POS_TICK_INTERVAL):
            angle = np.radians(90 - (bp / self.seq_len) * 360)
            r0, r1 = R - 0.04, R + 0.04
            ax.plot([np.cos(angle) * r0, np.cos(angle) * r1],
                    [np.sin(angle) * r0, np.sin(angle) * r1],
                    color='#aaa', linewidth=1)
            ax.text(np.cos(angle) * (r1 + 0.07), np.sin(angle) * (r1 + 0.07),
                    str(bp), fontsize=7, ha='center', va='center', color='#aaa')

    # --- Stacking constants --------------------------------------------
    MIN_ANG_GAP  : float = 0.10   # minimum angular gap (radians) before stacking
    STACK_STEP   : float = 0.26   # radial step per stack level
    STEM_COLOR   : str   = '#aaaaaa'  # colour of the leader line from tick to label

    def _compute_label_radii(self) -> list:
        """
        Purpose:
        --------
            Calculate the radius at which each enzyme label should be drawn. 
            - Labels that are angularly too close to an already-placed label at the same radial band are pushed outward by STACK_STEP until they fit
        Returns:
        --------
            Sorted list of tuples, one per cut event:
                (angle_rad, enzyme_name, position, colour, r_tick, r_label)
            - angle_rad : float
                Angular position on the circle in radians.
            - enzyme_name : str
                Name of the restriction enzyme.
            - position : int
                Absolute sequence position of the cut site.
            - colour : tuple
                RGB colour assigned to this enzyme.
            - r_tick : float
                Radius at which the tick mark ends.
            - r_label : float
                Radius at which the label text should be anchored, pushed
                outward as needed to avoid overlap with nearby labels.
        """
        R = self.RADIUS
        events = []
        for enzyme, data in self.cutters.items():
            col   = self.color_map[enzyme]
            multi = data[2] > 1
            for pos in data[3]:
                ang    = np.radians(90 - (pos / self.seq_len) * 360)
                r_tick = R + self.TICK_R + (self.MC_GAP if multi else 0)
                events.append((ang, enzyme, pos, col, r_tick))

        events.sort(key=lambda e: e[0])

        placed = []   # (angle, r_label) of already-positioned labels
        result = []
        for ang, enzyme, pos, col, r_tick in events:
            r_lbl = r_tick + self.LBL_GAP
            while True:
                conflict = any(
                    abs(p_r - r_lbl) < 0.05 and
                    min(abs(ang - p_a), 2*math.pi - abs(ang - p_a)) < self.MIN_ANG_GAP
                    for p_a, p_r in placed
                )
                if not conflict:
                    break
                r_lbl += self.STACK_STEP
            placed.append((ang, r_lbl))
            result.append((ang, enzyme, pos, col, r_tick, r_lbl))

        return result

    def _draw_cut_sites(self, ax) -> None:
        """
        Purpose:
        --------
            Draw a radial tick, a leader line, and a stacked label for every cut site.
            - Uses _compute_lable_radii() list to prevent label overlap.

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The axes object on which all cut-site annotations are drawn.
        """
        R       = self.RADIUS
        events  = self._compute_label_radii()

        for ang, enzyme, pos, col, r_tick, r_lbl in events:
            cos_a = np.cos(ang)
            sin_a = np.sin(ang)

            # Radial tick on the backbone
            ax.plot([cos_a * (R - self.TICK_R / 2), cos_a * r_tick],
                    [sin_a * (R - self.TICK_R / 2), sin_a * r_tick],
                    color=col, linewidth=2.5, solid_capstyle='round', zorder=4)

            # Leader line from tick tip to label anchor
            ax.plot([cos_a * r_tick, cos_a * r_lbl],
                    [sin_a * r_tick, sin_a * r_lbl],
                    color=self.STEM_COLOR, linewidth=0.8,
                    linestyle='--', alpha=0.6, zorder=3)

            # Label
            ha = 'left' if cos_a >= 0 else 'right'
            va = 'bottom' if sin_a >= 0 else 'top'
            ax.text(
                cos_a * r_lbl, sin_a * r_lbl,
                f"{enzyme}\n({pos})",
                fontsize=7.5, ha=ha, va=va,
                color=col, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.15', fc='white', ec=col, alpha=0.90, linewidth=0.7),
                zorder=5,
            )

    def _draw_centre_text(self, ax) -> None:
        """
        Purpose:
        --------
            Draw title, bp count, and cutter count at the centre of the circle.
        """
        ax.text(0,  0.10, self.title,
                ha='center', va='center', fontsize=14, fontweight='bold', color='#222')
        ax.text(0, -0.10, f"{self.seq_len:,} bp",
                ha='center', va='center', fontsize=11, color='#555')
        ax.text(0, -0.28, f"{len(self.cutters)} unique cutters",
                ha='center', va='center', fontsize=9, color='#888')

    def _draw_legend(self, ax) -> None:
        """
        Purpose:
        --------
            Draw the enzyme colour legend and a non-cutter list note.
        """
        handles = [
            mpatches.Patch(color=self.color_map[n], label=n)
            for n in self.cutters
        ]
        if handles:
            ax.legend(handles=handles, loc='lower left', fontsize=8, title='Cutting enzymes', title_fontsize=8.5, framealpha=0.9, ncol=1, bbox_to_anchor=(-0.22, 0))

        non_cutters = [k for k, v in self.results.items() if v[2] == 0]
        if non_cutters:
            ax.get_figure().text(
                0.5, 0.01,
                f"No cut: {', '.join(non_cutters)}",
                ha='center', va='bottom', fontsize=8, color='#999',
                wrap=True,
            )

    # --- Public entry point ---------------------------------------------

    def render(self) -> plt.Figure:
        """
        Purpose:
        --------
            Build and save the circular map.
        
        Returns:
        --------
            fig : matplotlib.figure.Figure
                The fully rendered figure object. Written to disk at
                self.OUTPUT_PATH.
        """
        sys.stdout.write('\n\t2. Generating circular plasmid map...\n')
        fig, ax = self._make_figure()
        self._draw_backbone(ax)
        self._draw_position_ticks(ax)
        self._draw_cut_sites(ax)
        self._draw_centre_text(ax)
        self._draw_legend(ax)
        plt.tight_layout()
        if self.OUTPUT_PATH:
            fig.savefig(self.OUTPUT_PATH, dpi=self.DPI, bbox_inches='tight')
            sys.stdout.write(f'\t--> Circular map saved to {HIGHLIGHT}{self.OUTPUT_PATH}{RESET}\n')
        return fig


# =======================================================
# LinearMap  (single-stranded sequence viewer)
# =======================================================

class LinearMap(_EnzymeDataMixin):
    """
    Single-stranded sequence viewer.  Wraps the sequence into lines of
    CHARS_PER_LINE bases, colours each base by identity, highlights recognition
    motifs, and marks top-strand cut sites with a vertical tick in the inter-base
    gap plus a labelled arrow.

    Parameters
    ----------
    results : dict
        { enzyme_name: [recognition_seq, cut_notation, cut_count, [site_positions]] }
    plasmid_sequence : str
        Full nucleotide sequence (top strand, 5'→3').
    title : str
        Header text.
    figwidth : float
        Figure width in inches.
    chars_per_line : int
        Bases per wrapped line.
    dpi : int
        Output resolution.
    output_path : str
        Destination file.
    """

    # --- Geometry constants --------------------------------------------
    SEQ_LINE_H : float = 0.62   # inches per wrapped line
    LEFT_MAR   : float = 0.065  # x-fraction reserved for position labels
    HEADER_H   : float = 1.4    # fixed height of the header panel in inches
    ARROW_BASE : float = 0.3    # minimum arrow extension above/below tick
    LEVEL_STEP : float = 0.8   # extra arrow height per horizontal stack level
    TICK_HALF  : float = 0.44   # half-height of the cut-site tick mark
    MIN_X_GAP  : float = 0.08   # minimum x-distance between labels before stacking

    def __init__(
        self,
        results: dict,
        plasmid_sequence: str,
        title: str = "Plasmid",
        figwidth: float = 14.0,
        chars_per_line: int = 80,
        dpi: int = 150,
        output_path: str = "results/ss_linear_map.png",
    ) -> None:
        self.results        = results
        self.seq            = plasmid_sequence.upper().replace(' ', '').replace('\n', '')
        self.seq_len        = len(self.seq)
        self.title          = title
        self.FIGWIDTH       = figwidth
        self.CHARS_PER_LINE = chars_per_line
        self.DPI            = dpi
        self.OUTPUT_PATH    = output_path

        self.enzyme_names  = list(results.keys())
        self.color_map     = self._build_color_map()
        self.char_enzyme   = self._build_char_enzyme_map()
        self.cut_positions = self._build_top_cut_positions()

    # --- Figure setup -------------------------------------------------

    def _cuts_per_line(self) -> dict:
        """
        Purpose:
        --------
            Helper function to determine how many sequence sections contain enzyme cuts. 
            - used by _make_figure() to determine how much extra vertical space to add to the figure for cut arrows.
            - parse self.cut_positions lookup table from _build_top_cut_positions()  
                - {top cut index: enzyme name list}
            - Determine max enzyme stack depth out of all lines.
                - example:
                    position 283 → ['EcoRI']           line 3, stack depth 1
                    position 262 → ['BamHI']           line 3, stack depth 1
                    position 35  → ['HaeIII', 'SmaI']  line 0, stack depth 2
                    - returns: {3: 1, 0: 2}
                    
        Returns:
        --------
            {line_idx: max_stack_depth} for lines containing cut sites.
                - line_idx : int
                    Zero-based index of the wrapped sequence line.
                - max_stack_depth : int
                    The highest horizontal stack level assigned to any cut site
                    on that line, used to calculate extra vertical clearance.
        """
        char_w    = (1.0 - self.LEFT_MAR - 0.005) / self.CHARS_PER_LINE
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        result: dict[int, int] = {}
        for li in range(num_lines):
            start = li * self.CHARS_PER_LINE
            end   = min(start + self.CHARS_PER_LINE, self.seq_len)
            levels = self._assign_cut_levels(self.cut_positions, start, end, char_w)
            if levels:
                result[li] = max(levels.values())
        return result

    def _make_figure(self) -> tuple:
        """
        Purpose:
        --------
            Create the figure with a fixed-height header panel and a
            dynamically-sized sequence panel. The sequence panel height scales
            with the number of wrapped lines plus the vertical clearance required
            for cut-site arrow annotations.

        Returns:
        --------
            fig : matplotlib.figure.Figure
                The top-level figure object with a dark (#1e2228) background.
            ax_hdr : matplotlib.axes.Axes
                The fixed-height header axes for the title, bp count, and legend.
            ax_seq : matplotlib.axes.Axes
                The variable-height sequence axes for all wrapped base rows and
                cut-site annotations.
        """

        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        # Each level adds LEVEL_STEP=0.35 data-units of arrow height.
        # Convert to figure inches using the seq panel's data-unit-to-inch ratio.

        extra     = sum((d + 1) * self.LEVEL_STEP * self.SEQ_LINE_H
                        for d in self._cuts_per_line().values())
        seq_h     = num_lines * self.SEQ_LINE_H + extra
        fig_h     = self.HEADER_H + seq_h

        fig = plt.figure(figsize=(self.FIGWIDTH, fig_h))
        fig.patch.set_facecolor('#1e2228')

        gs = fig.add_gridspec(2, 1,
                            height_ratios=[self.HEADER_H, seq_h],
                            hspace=0.04)
        ax_hdr = fig.add_subplot(gs[0])
        ax_seq = fig.add_subplot(gs[1])
        for ax in (ax_hdr, ax_seq):
            ax.set_facecolor('#1e2228')
        return fig, ax_hdr, ax_seq

    # --- Drawing helpers ---------------------------------------------

    def _draw_header(self, ax) -> None:
        """
        Purpose:
        --------
            Draw the title, bp count, and enzyme colour legend in the dedicated
            header panel.

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The header axes object. Expected to span the full figure width
                with normalised coordinates [0, 1] x [0, 1].
        """
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        ax.text(0.0, 0.95, self.title,
                ha='left', va='top', fontsize=14, fontweight='700',
                color='#00d4ff', fontfamily='monospace',
                transform=ax.transAxes)
        ax.text(1.0, 0.95, f"{self.seq_len} bp",
                ha='center', va='top', fontsize=12,
                color='#00d4ff', fontfamily='monospace',
                transform=ax.transAxes)
        handles = self._build_legend_handles(self.color_map)
        if handles:
            leg = ax.legend(
                handles=handles,
                loc='upper left', bbox_to_anchor=(0.0, 0.80),
                bbox_transform=ax.transAxes,
                ncol=min(4, len(handles)),
                fontsize=9, frameon=True, framealpha=0.12,
                edgecolor='#2a4a6b', facecolor='#252b33', labelcolor='#cdd9e5',
                handlelength=1.0, handleheight=0.8, columnspacing=1.2,
            )
            for t in leg.get_texts():
                t.set_fontfamily('monospace')

    def _assign_cut_levels(self, cut_positions: dict, start: int, end: int, char_w: float, min_x_gap: float = None) -> dict:
        """
        Purpose:
        --------
        Assign a vertical stack level to each cut site on one wrapped line so
        that labels never collide horizontally.

        Works left-to-right along the line. Each cut site starts at level 0.
        If its x position is within min_x_gap of any already-placed site at
        the same level, it is bumped to the next level up.  Levels translate
        directly into arrow heights in _draw_cut_annotation.

        Arguments:
        ----------
        cut_positions : dict
            The full cut_positions dict (abs_idx -> [enzyme names]).
        start : int
            Absolute sequence index of the first base on this wrapped line.
        end : int
            Absolute sequence index one past the last base on this wrapped line.
        char_w : float
            Width of one character in normalised axis units.
        min_x_gap : float, optional
            Minimum x-distance between two labels at the same level before
            the newer one is bumped upward. Defaults to self.MIN_X_GAP.

        Returns:
        --------
            dict[int, int]
                Mapping of absolute sequence index to stack level (0 = lowest /
                closest to the sequence row). Only indices within [start, end)
                that appear in cut_positions are included.
        """

        if min_x_gap is None:
            min_x_gap = self.MIN_X_GAP
        sites = sorted(
            (idx for idx in cut_positions if start <= idx < end),
            key=lambda idx: idx
        )
        levels: dict[int, int] = {}
        level_xs: list[list[float]] = []  # all x positions per level, not just last

        for idx in sites:
            x = self.LEFT_MAR + (idx - start) * char_w

            placed = False
            for lv, xs in enumerate(level_xs):
                if all(abs(x - px) > min_x_gap for px in xs):  # check against ALL at this level
                    levels[idx] = lv
                    xs.append(x)
                    placed = True
                    break
            if not placed:
                levels[idx] = len(level_xs)
                level_xs.append([x])

        return levels

    def _line_y_positions(self) -> list:
        """
        Purpose:
        --------
            Compute the y data-coordinate for each wrapped sequence line, accounting
            for the vertical clearance consumed by arrows on the preceding line.

            Each line's y is offset downward by the maximum arrow height of the line
            above it, so annotations never intrude into the sequence row below.

        Returns:
        --------
            list of float
                Y data-coordinates for each wrapped line, in top-to-bottom order
                (values decrease as line index increases). The list has one entry
                per wrapped line (length = ceil(seq_len / CHARS_PER_LINE)).
        """

        char_w    = (1.0 - self.LEFT_MAR - 0.005) / self.CHARS_PER_LINE
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        cpl       = self._cuts_per_line()

        y_positions = []
        y = -0.5
        for li in range(num_lines):
            # Apply THIS line's clearance before placing it, not after
            max_level = cpl.get(li, 0)
            clearance = self.TICK_HALF + self.ARROW_BASE + max_level * self.LEVEL_STEP + 0.10
            y -= clearance
            y_positions.append(y)
            y -= 1.0  # base line spacing
        return y_positions


    def _setup_seq_axes(self, ax, num_lines: int) -> None:
        """
        Purpose:
        --------
            Configure the coordinate space for the sequence panel by setting
            axis limits that accommodate all wrapped lines and their annotations.

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The sequence axes object to configure.
            num_lines : int
                Total number of wrapped sequence lines to be rendered.
        """
        y_positions = self._line_y_positions()
        y_min = y_positions[-1] - 0.6
        ax.set_xlim(0, 1)
        ax.set_ylim(y_min, 0.5)
        ax.axis('off')

    def _draw_all_lines(self, ax) -> None:
        """
        Purpose:
        --------
            Render every wrapped sequence line by iterating over all line indices
            and delegating each to _draw_sequence_line().

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The sequence axes object on which all lines are drawn.
        """

        char_w      = (1.0 - self.LEFT_MAR - 0.005) / self.CHARS_PER_LINE
        y_positions = self._line_y_positions()
        num_lines   = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        for li in range(num_lines):
            self._draw_sequence_line(ax, li, char_w, y_positions[li])

    def _draw_sequence_line(self, ax, line_idx: int, char_w: float, y: float = None) -> None:
        """
        Purpose:
        --------
            Render a single wrapped sequence row: the position label on the left,
            each nucleotide character with optional motif highlight, and any
            cut-site tick marks and arrows that fall within this line's range.

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The sequence axes object on which this line is drawn.
            line_idx : int
                Zero-based index of the wrapped line to render.
            char_w : float
                Width of one character in normalised axis units.
            y : float, optional
                Y data-coordinate at which to centre the sequence row. Defaults
                to -(line_idx + 0.5) if not provided.
        """
        
        if y is None:
            y = -(line_idx + 0.5)
        start = line_idx * self.CHARS_PER_LINE
        end   = min(start + self.CHARS_PER_LINE, self.seq_len)

        ax.text(self.LEFT_MAR - 0.004, y, f"{start + 1:>6}",
                ha='right', va='center', fontsize=7,
                color='#3a6a9b', fontfamily='monospace')

        cut_levels = self._assign_cut_levels(self.cut_positions, start, end, char_w)

        for ci, base in enumerate(self.seq[start:end]):
            abs_idx = start + ci
            x_left  = self.LEFT_MAR + ci * char_w
            x_mid   = x_left + char_w / 2
            self._draw_base(ax, base, abs_idx, x_left, x_mid, y, char_w)
            if abs_idx in self.cut_positions:
                self._draw_cut_annotation(ax, abs_idx, x_left, y, char_w,
                                        self.cut_positions, above=True,
                                        level=cut_levels.get(abs_idx, 0))


    def _draw_base(self, ax, base: str, abs_idx: int, x_left: float, x_mid: float, y: float, char_w: float) -> None:
        """
        Purpose:
        --------
            Draw one nucleotide character with an optional coloured highlight
            background when the base falls within a recognised enzyme motif.
            Bases inside a motif are rendered white and bold; bases outside use
            their identity colour from BASE_COLORS.

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The sequence axes on which the base is drawn.
            base : str
                Single nucleotide character ('A', 'T', 'G', or 'C').
            abs_idx : int
                Absolute sequence index of this base, used to look up motif
                membership in self.char_enzyme.
            x_left : float
                Left edge x-coordinate of this character's cell in normalised
                axis units.
            x_mid : float
                Centre x-coordinate of this character's cell.
            y : float
                Y data-coordinate of the sequence row centre.
            char_w : float
                Width of one character cell in normalised axis units.
        """
        hit = self.char_enzyme.get(abs_idx)
        if hit:
            ax.add_patch(mpatches.FancyBboxPatch(
                (x_left, y - 0.40), char_w, 0.80,
                boxstyle='square,pad=0.0', linewidth=0,
                facecolor=self.color_map[hit], alpha=0.35, zorder=2,
            ))
            txt_col, bold = '#ffffff', True
        else:
            txt_col, bold = BASE_COLORS.get(base, '#cdd9e5'), False

        ax.text(x_mid, y, base, ha='center', va='center', fontsize=7.5,
                color=txt_col, fontfamily='monospace',
                fontweight='bold' if bold else 'normal', zorder=3)


    def _draw_cut_annotation(self, ax, abs_idx: int, x_left: float, y: float, char_w: float, cut_positions: dict, above: bool = True, level: int = 0) -> None:
        """
        Purpose:
        --------
            Draw the cut-site annotation for one inter-base gap: a vertical tick
            spanning the sequence row and an upward or downward arrow with an
            enzyme name and position label.

            The tick is centred on x_left (the left edge of the cut base).
            Arrows point away from the sequence row: upward for top-strand cuts
            (above=True), downward for bottom-strand cuts (above=False).
            Multiple enzymes at the same position are stacked vertically.
            The level parameter adds extra height so nearby labels on the same
            line do not overlap horizontally.

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The sequence axes on which the annotation is drawn.
            abs_idx : int
                Absolute sequence index of the cut position; used as a key into
                cut_positions to retrieve the list of enzyme names.
            x_left : float
                X-coordinate of the inter-base gap (left edge of the cut base)
                in normalised axis units.
            y : float
                Y data-coordinate of the sequence strand row.
            char_w : float
                Width of one character cell, used to offset the label text
                slightly to the right of the tick.
            cut_positions : dict
                Lookup table {abs_idx: [enzyme_name, ...]} for the strand being
                annotated (top or bottom).
            above : bool, optional
                If True (default), arrows point upward (top strand). If False,
                arrows point downward (bottom strand).
            level : int, optional
                Horizontal stack level from _assign_cut_levels(). Higher levels
                produce taller arrows to clear neighbouring labels. Default is 0.
        """

        direction  = 1 if above else -1

        for stack_i, enz_name in enumerate(cut_positions[abs_idx]):
            col      = self.color_map[enz_name]
            tick_top = y + self.TICK_HALF
            tick_bot = y - self.TICK_HALF
            tip_y    = tick_top if above else tick_bot
            origin_y = tip_y + direction * (0.60 + stack_i * 0.30 + level * self.LEVEL_STEP)

            ax.plot([x_left, x_left], [tick_bot, tick_top],
                    color=col, lw=2.0, zorder=6, solid_capstyle='round')
            ax.annotate('',
                        xy=(x_left, tip_y), xytext=(x_left, origin_y),
                        arrowprops=dict(arrowstyle='->', color=col, lw=1.1),
                        zorder=5)
            ax.text(x_left + char_w * 0.3, origin_y,
                    f"{enz_name} ↓{abs_idx}",
                    ha='left', va='center', fontsize=6.5,
                    color=col, fontfamily='monospace',
                    fontweight='600', zorder=5)

    # --- Public entry point --------------------------------------------

    def render(self) -> plt.Figure:
        """
        Purpose:
        --------
            Build, display, and save the single-stranded linear map.

        Returns:
        --------
            fig : matplotlib.figure.Figure
                The fully rendered figure object. Also written to disk at
                self.OUTPUT_PATH if that attribute is set.
        """
        sys.stdout.write('\n\t3. Generating single-stranded linear plasmid map...\n')
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        fig, ax_hdr, ax_seq = self._make_figure()
        self._draw_header(ax_hdr)
        self._setup_seq_axes(ax_seq, num_lines)
        self._draw_all_lines(ax_seq)
        fig.subplots_adjust(left=0.03, right=0.97, top=0.97, bottom=0.03, hspace=0.04)
        if self.OUTPUT_PATH:
            fig.savefig(self.OUTPUT_PATH, dpi=self.DPI, bbox_inches='tight',
                        facecolor=fig.get_facecolor())
            sys.stdout.write(f'\t--> Linear map saved to {HIGHLIGHT}{self.OUTPUT_PATH}{RESET}\n')
        return fig


# =======================================================
# DoubleStrandedMap
# =======================================================

class DoubleStrandedMap(LinearMap):
    """
    Double-stranded sequence viewer.  Inherits all single-strand drawing from
    LinearMap and adds a bottom (complement, 3'→5') strand beneath each top
    strand line.

    Layout per wrapped block
    ------------------------
        [annotation arrows — top strand cuts — pointing UP]
        5'  A T G C A T …  3'   ← top strand, coloured by base
            | | | | | |         ← pairing tick marks
        3'  T A C G T A …  5'   ← bottom strand, coloured by base (complement)
        [annotation arrows — bottom strand cuts — pointing DOWN]

    Top-strand cut positions use the ^ offset from the notation.
    Bottom-strand cut positions use the _ offset from the notation.
    Both arrows are drawn in the inter-base gap (x_left of the cut index).
    """

    # Geometry for within-block layout (strand rows + pairing marks + arrows)
    BOT_STRAND_H  : float = 0.55   # inches of extra figure height per block
    PAIR_TICK_H   : float = 0.15   # half-height of base-pairing tick marks
    TOP_ARROW_H   : float = 0.44   # data-units above top strand for upward arrows
    BOT_ARROW_H   : float = 0.44   # data-units below bottom strand for downward arrows
    STRAND_GAP    : float = 0.70   # data-unit gap between top and bottom strand rows
    BLOCK_GAP     : float = 1.20   # extra data-units of whitespace *between* blocks
    LEVEL_STEP    : float = 0.55   # override LinearMap's value just for ds map

    def __init__(
        self,
        results: dict,
        plasmid_sequence: str,
        title: str = "Plasmid",
        figwidth: float = 14.0,
        chars_per_line: int = 80,
        dpi: int = 150,
        output_path: str = "results/ds_linear_map.png",
    ) -> None:
        super().__init__(results, plasmid_sequence, title, figwidth, chars_per_line, dpi, output_path)

        # Complement strand (3'→5', so same left-to-right order as top strand)
        self.comp_seq       = ''.join(COMPLEMENT.get(b, 'N') for b in self.seq)
        # Bottom-strand cut positions (keyed by the same abs index as top strand)
        self.bot_cut_positions = self._build_bot_cut_positions()

    # --- Override geometry ---------------------------------------------

    def _seq_line_height(self) -> float:
        """
        Purpose:
        --------
            Return the total figure height in inches consumed by one
            double-stranded block, combining the top-strand row height
            with the additional space required for the bottom strand.

        Returns:
        --------
            float
                Combined height in inches for one double-stranded block
                (SEQ_LINE_H + BOT_STRAND_H).
        """
        return self.SEQ_LINE_H + self.BOT_STRAND_H

    def _cuts_per_line_bot(self) -> dict:
        """
        Purpose:
        --------
            Compute the maximum horizontal stack depth per wrapped line for
            bottom-strand cut-site annotations.

        Returns:
        --------
            {line_idx: max_stack_depth}
                Dictionary mapping zero-based line index to the highest stack
                level assigned to any bottom-strand cut site on that line.
                Lines with no bottom-strand cuts are omitted.
        """
        char_w    = (1.0 - self.LEFT_MAR - 0.005) / self.CHARS_PER_LINE
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        result: dict[int, int] = {}
        for li in range(num_lines):
            start = li * self.CHARS_PER_LINE
            end   = min(start + self.CHARS_PER_LINE, self.seq_len)
            levels = self._assign_cut_levels(self.bot_cut_positions, start, end, char_w)
            if levels:
                result[li] = max(levels.values())
        return result

    def _block_height(self, line_idx: int, top_cpl: dict, bot_cpl: dict) -> float:
        """
        Purpose:
        --------
            Calculate the total data-unit height of one double-stranded block,
            including the vertical clearance needed above the top strand for
            upward arrows and below the bottom strand for downward arrows.

        Arguments:
        ----------
            line_idx : int
                Zero-based index of the wrapped line whose block height is needed.
            top_cpl : dict
                {line_idx: max_stack_depth} for top-strand cuts, as returned
                by _cuts_per_line().
            bot_cpl : dict
                {line_idx: max_stack_depth} for bottom-strand cuts, as returned
                by _cuts_per_line_bot().

        Returns:
        --------
            float
                Total height in data units for this block:
                top_clearance + STRAND_GAP + bottom_clearance.
                Clearance for each strand is TICK_HALF + ARROW_BASE +
                (stack_depth * LEVEL_STEP) + 0.10, or 0.10 if no cuts exist
                on that strand for this line.
        """

        top_levels = top_cpl.get(line_idx, 0)
        bot_levels = bot_cpl.get(line_idx, 0)
        top_clear  = (self.TICK_HALF + self.ARROW_BASE + top_levels * self.LEVEL_STEP + 0.10) if line_idx in top_cpl else 0.10
        bot_clear  = (self.TICK_HALF + self.ARROW_BASE + bot_levels * self.LEVEL_STEP + 0.10) if line_idx in bot_cpl else 0.10
        return top_clear + self.STRAND_GAP + bot_clear

    def _make_figure(self) -> tuple:
        """
        Purpose:
        --------
            Create the figure with a fixed-height header panel and a
            dynamically-sized sequence panel sized for double-stranded blocks.
            Total sequence panel height is the sum of all per-block heights
            plus inter-block gap spacing.

        Returns:
        --------
            fig : matplotlib.figure.Figure
                The top-level figure object with a dark (#1e2228) background.
            ax_hdr : matplotlib.axes.Axes
                The fixed-height header axes for the title, bp count, and legend.
            ax_seq : matplotlib.axes.Axes
                The variable-height sequence axes for all double-stranded blocks
                and their cut-site annotations.
        """

        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        top_cpl   = self._cuts_per_line()
        bot_cpl   = self._cuts_per_line_bot()
        seq_h     = sum(self._block_height(li, top_cpl, bot_cpl) + self.BLOCK_GAP
                        for li in range(num_lines))
        fig_h     = self.HEADER_H + seq_h

        fig = plt.figure(figsize=(self.FIGWIDTH, fig_h))
        fig.patch.set_facecolor('#1e2228')
        gs = fig.add_gridspec(2, 1, height_ratios=[self.HEADER_H, seq_h], hspace=0.04)
        ax_hdr = fig.add_subplot(gs[0])
        ax_seq = fig.add_subplot(gs[1])
        for ax in (ax_hdr, ax_seq):
            ax.set_facecolor('#1e2228')
        return fig, ax_hdr, ax_seq

    def _setup_seq_axes(self, ax, num_lines: int) -> None:
        """
        Purpose:
        --------
            Configure the coordinate space for the sequence panel to accommodate
            all double-stranded blocks and their annotations.
            - Modifies ax in place by setting xlim to [0, 1], ylim to [y_min, 0.5] where y_min is calculated from the cumulative height of all blocks and inter-block gaps
            - turns off axis decorations.

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The sequence axes object to configure.
            num_lines : int
                Total number of wrapped double-stranded blocks to be rendered.

        Returns:
        --------
            None. 
        """
        top_cpl = self._cuts_per_line()
        bot_cpl = self._cuts_per_line_bot()
        y_min   = -sum(self._block_height(li, top_cpl, bot_cpl) + self.BLOCK_GAP
                    for li in range(num_lines)) - 0.5
        ax.set_xlim(0, 1)
        ax.set_ylim(y_min, 0.5)
        ax.axis('off')

    # --- Override line drawing -----------------------------------------

    def _draw_all_lines(self, ax) -> None:
        """
        Purpose:
        --------
            Render every wrapped double-stranded block by iterating over all
            line indices, computing per-block y positions, and delegating
            each block to _draw_ds_block().

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The sequence axes object on which all blocks are drawn.
        """
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        char_w    = (1.0 - self.LEFT_MAR - 0.005) / self.CHARS_PER_LINE
        top_cpl   = self._cuts_per_line()
        bot_cpl   = self._cuts_per_line_bot()

        y = -0.5
        for li in range(num_lines):
            top_levels = top_cpl.get(li, 0)
            bot_levels = bot_cpl.get(li, 0)
            
            # Only add clearance if this line actually has cuts
            top_clear = (self.TICK_HALF + self.ARROW_BASE + top_levels * self.LEVEL_STEP + 0.10) if top_levels > 0 or li in top_cpl else 0.10
            bot_clear = (self.TICK_HALF + self.ARROW_BASE + bot_levels * self.LEVEL_STEP + 0.10) if bot_levels > 0 or li in bot_cpl else 0.10

            y -= top_clear
            self._draw_ds_block(ax, li, char_w, y_top=y)
            y -= self.STRAND_GAP + bot_clear + self.BLOCK_GAP

    def _draw_ds_block(self, ax, line_idx: int, char_w: float, y_top: float = None) -> None:
        """
        Purpose:
        --------
            Draw one double-stranded block consisting of a top strand row and a
            complementary bottom strand row, connected by base-pairing tick marks,
            with strand polarity labels (5'/3') and cut-site annotations for both
            strands.
            Modifies ax in place by adding:
                - Position label and 5'/3' polarity labels on both strands.
                - Per-base character glyphs with optional motif highlights for both the top (sense) strand and the bottom (complement) strand.
                - Base-pairing tick marks between each paired nucleotide.
                - Upward cut-site annotations for top-strand cut positions.
                - Downward cut-site annotations for bottom-strand cut positions.

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The sequence axes object on which the block is drawn.
            line_idx : int
                Zero-based index of the wrapped line to render.
            char_w : float
                Width of one character cell in normalised axis units.
            y_top : float, optional
                Y data-coordinate for the centre of the top-strand row. If not
                provided, it is calculated from line_idx using the standard
                block height and TOP_ARROW_H offset.
        """
        if y_top is None:
            line_h = self._seq_line_height()
            y_top  = -(line_idx * (line_h + self.BLOCK_GAP) + self.TOP_ARROW_H)
        y_bot   = y_top - self.STRAND_GAP

        start = line_idx * self.CHARS_PER_LINE
        end   = min(start + self.CHARS_PER_LINE, self.seq_len)

        # Compute stack levels for both strands
        top_levels = self._assign_cut_levels(self.cut_positions, start, end, char_w)
        bot_levels = self._assign_cut_levels(self.bot_cut_positions, start, end, char_w)

        y_mid = (y_top + y_bot) / 2
        ax.text(self.LEFT_MAR - 0.004, y_mid, f"{start + 1:>6}",
                ha='right', va='center', fontsize=7,
                color='#3a6a9b', fontfamily='monospace')

        ax.text(self.LEFT_MAR - 0.004, y_top, "5'", ha='right', va='bottom', fontsize=6, color='#7a96aa', fontfamily='monospace')
        ax.text(self.LEFT_MAR - 0.004, y_bot, "3'", ha='right', va='top',    fontsize=6, color='#7a96aa', fontfamily='monospace')
        ax.text(1.0,                   y_top, "3'", ha='left',  va='bottom', fontsize=6, color='#556677', fontfamily='monospace')
        ax.text(1.0,                   y_bot, "5'", ha='left',  va='top',    fontsize=6, color='#556677', fontfamily='monospace')

        for ci in range(end - start):
            abs_idx  = start + ci
            top_base = self.seq[abs_idx]
            bot_base = self.comp_seq[abs_idx]
            x_left   = self.LEFT_MAR + ci * char_w
            x_mid    = x_left + char_w / 2

            self._draw_base(ax, top_base, abs_idx, x_left, x_mid, y_top, char_w)
            self._draw_base(ax, bot_base, abs_idx, x_left, x_mid, y_bot, char_w)
            self._draw_pair_tick(ax, x_mid, y_top, y_bot)

            if abs_idx in self.cut_positions:
                self._draw_cut_annotation(ax, abs_idx, x_left, y_top, char_w,
                                        self.cut_positions, above=True,
                                        level=top_levels.get(abs_idx, 0))
            if abs_idx in self.bot_cut_positions:
                self._draw_cut_annotation(ax, abs_idx, x_left, y_bot, char_w,
                                        self.bot_cut_positions, above=False,
                                        level=bot_levels.get(abs_idx, 0))

        for ci in range(end - start):
            self._draw_pair_tick(ax, self.LEFT_MAR + ci * char_w + char_w / 2, y_top, y_bot)

    def _draw_pair_tick(self, ax, x_mid: float, y_top: float, y_bot: float) -> None:
        """
        Purpose:
        --------
            Draw a short vertical line between a paired base pair to indicate
            Watson-Crick hydrogen bonding between the top and bottom strands.

        Arguments:
        ----------
            ax : matplotlib.axes.Axes
                The sequence axes on which the tick is drawn.
            x_mid : float
                X data-coordinate at the horizontal centre of the base pair column.
            y_top : float
                Y data-coordinate of the top-strand row centre.
            y_bot : float
                Y data-coordinate of the bottom-strand row centre.

        """
        gap    = y_bot - y_top          # negative value
        centre = y_top + gap / 2
        ax.plot([x_mid, x_mid],
                [centre - self.PAIR_TICK_H, centre + self.PAIR_TICK_H],
                color='#334455', lw=0.6, zorder=1)

    # --- Public entry point --------------------------------------------

    def render(self) -> plt.Figure:
        """
        Purpose:
        --------
            Build, display, and save the double-stranded linear map.

        Returns:
        --------
            fig : matplotlib.figure.Figure
                The fully rendered figure object. Also written to disk at
                self.OUTPUT_PATH
        """
        sys.stdout.write('\n\t3. Generating double-stranded linear plasmid map...\n')
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        fig, ax_hdr, ax_seq = self._make_figure()
        self._draw_header(ax_hdr)
        self._setup_seq_axes(ax_seq, num_lines)
        self._draw_all_lines(ax_seq)
        fig.subplots_adjust(left=0.03, right=0.97, top=0.97, bottom=0.03, hspace=0.04)
        if self.OUTPUT_PATH:
            fig.savefig(self.OUTPUT_PATH, dpi=self.DPI, bbox_inches='tight',
                        facecolor=fig.get_facecolor())
            sys.stdout.write(f'\t--> Double-stranded map saved to {HIGHLIGHT}{self.OUTPUT_PATH}{RESET}\n')
        return fig

class PlasmidMap:
    """
    Convenience wrapper that exposes circular and linear rendering through a
    single object.  All heavy lifting is delegated to CircularMap, LinearMap,
    and DoubleStrandedMap.
    """

    def __init__(self, results: dict, plasmid_sequence: str, title: str = "Plasmid name") -> None:
        self.results          = results
        self.plasmid_sequence = plasmid_sequence
        self.title            = title

    def annotate_circular(self, figsize: tuple = (10, 10), dpi: int = 150, output_path: str = "results/circular_map.png") -> plt.Figure:
        """
        Purpose:
        --------
            Render the circular plasmid map.
        """
        short_title = self.title.split(",")[0]
        return CircularMap(
            self.results, self.plasmid_sequence, short_title, 
            figsize=figsize, 
            dpi=dpi, 
            output_path=output_path
        ).render()

    def annotate_linear(self, figwidth: float = 14.0, chars_per_line: int = 80, dpi: int = 150, output_path: str = "results/ss_linear_map.png") -> plt.Figure:
        """
        Purpose:
        --------
        Render a single-stranded linear sequence map.
        """
        return LinearMap(
            self.results, self.plasmid_sequence, self.title,
            figwidth=figwidth,
            chars_per_line=chars_per_line,
            dpi=dpi,
            output_path=output_path,
        ).render()

    def annotate_double_stranded(self, figwidth: float = 14.0, chars_per_line: int = 80, dpi: int = 150, output_path: str = "results/ds_linear_map.png") -> plt.Figure:
        """
        Purpose:
        --------
        Render a double-stranded linear sequence map.
        """
        return DoubleStrandedMap(
            self.results, self.plasmid_sequence, self.title,
            figwidth=figwidth,
            chars_per_line=chars_per_line,
            dpi=dpi,
            output_path=output_path,
        ).render()