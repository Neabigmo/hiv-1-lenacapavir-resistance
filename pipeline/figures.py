"""
Visualization Module for HIV Capsid Surveillance Analysis

Creates publication-quality figures including:
- World map with country coloring
- Bar charts for regional data
- Heatmaps for site coverage
- Radar charts for regional profiles
- All figures at 300+ DPI

Author: HIV Capsid Surveillance Analysis Pipeline
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict

import numpy as np
import pandas as pd

# Matplotlib configuration for publication quality
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D

# Optional: geopandas for map visualization
try:
    import geopandas as gpd
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False

from . import config

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class FigureStyle:
    """Manages consistent styling for all figures."""

    def __init__(self):
        self.dpi = config.PLOT_SETTINGS.get('dpi', 300)
        self.fig_size = config.PLOT_SETTINGS.get('fig_size', (10, 8))
        self.font_size = config.PLOT_SETTINGS.get('font_size', 12)
        self.title_size = config.PLOT_SETTINGS.get('title_font_size', 14)

        # Apply global style
        plt.rcParams.update({
            'font.family': 'sans-serif',
            'font.size': self.font_size,
            'axes.titlesize': self.title_size,
            'axes.labelsize': self.font_size,
            'xtick.labelsize': self.font_size - 1,
            'ytick.labelsize': self.font_size - 1,
            'legend.fontsize': self.font_size - 1,
            'figure.dpi': self.dpi,
            'savefig.dpi': self.dpi,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1,
        })

    def get_color_palette(self) -> Dict[str, str]:
        """Get region color palette."""
        return config.PLOT_SETTINGS.get('color_palette', {
            'AFRO': '#E64A19',
            'AMRO': '#1976D2',
            'SEARO': '#388E3C',
            'WPRO': '#7B1FA2',
            'EURO': '#00796B',
            'EMRO': '#C2185B',
        })


class WorldMapPlotter:
    """Creates world map visualizations of surveillance readiness."""

    def __init__(self, style: FigureStyle):
        self.style = style

    def plot_cir_sri_map(self, data: pd.DataFrame, output_path: str,
                        title: str = "CIR-SRI Scores by Country"):
        """
        Create a world map colored by CIR-SRI scores.

        Args:
            data: DataFrame with 'country' and 'cir_sri_score' columns
            output_path: Output file path
            title: Plot title
        """
        if not HAS_GEOPANDAS:
            logger.warning("geopandas not available, skipping world map")
            return

        fig, ax = plt.subplots(1, 1, figsize=(16, 10))

        # Load world shapefile (use natural earth data)
        try:
            world = gpd.read_file(
                gpd.datasets.get_path('naturalearth_lowres')
            )
        except Exception as e:
            logger.error(f"Could not load world data: {e}")
            return

        # Merge data with world geometries
        # ISO 3-letter codes need to be mapped to 2-letter
        iso3_to_iso2 = {
            'ZAF': 'ZA', 'KEN': 'KE', 'BWA': 'BW', 'MWI': 'MW', 'ZWE': 'ZW',
            'UGA': 'UG', 'ZMB': 'ZM', 'TZA': 'TZ', 'NGA': 'NG', 'ETH': 'ET',
            'USA': 'US', 'BRA': 'BR', 'PER': 'PE', 'CAN': 'CA', 'ARG': 'AR',
            'MEX': 'MX', 'COL': 'CO', 'THA': 'TH', 'IND': 'IN', 'IDN': 'ID',
            'VNM': 'VN', 'MMR': 'MM', 'CHN': 'CN', 'AUS': 'AU', 'JPN': 'JP',
            'KOR': 'KR', 'PHL': 'PH', 'GBR': 'GB', 'DEU': 'DE', 'FRA': 'FR',
            'BEL': 'BE', 'ESP': 'ES', 'ITA': 'IT', 'NLD': 'NL', 'RUS': 'RU',
            'EGY': 'EG', 'IRN': 'IR', 'PAK': 'PK', 'SDN': 'SD',
        }

        # Create country code mapping
        world['iso2'] = world['iso_a3'].map(iso3_to_iso2)

        # Merge with our data
        country_scores = data.set_index('country')['cir_sri_score'].to_dict()
        world['cir_sri'] = world['iso2'].map(country_scores)

        # Plot
        cmap = LinearSegmentedColormap.from_list(
            'cir_sri_cmap', ['#d73027', '#fee08b', '#1a9850']
        )

        world.boundary.plot(ax=ax, linewidth=0.5, color='gray')
        world.plot(
            column='cir_sri',
            ax=ax,
            cmap=cmap,
            missing_kwds={'color': 'lightgray', 'edgecolor': 'gray', 'linewidth': 0.5},
            legend=True,
            legend_kwds={
                'label': 'CIR-SRI Score (0-10)',
                'orientation': 'horizontal',
                'shrink': 0.6,
                'pad': 0.02,
            }
        )

        ax.set_title(title, fontsize=self.style.title_size + 2, fontweight='bold')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_xlim(-180, 180)
        ax.set_ylim(-60, 85)

        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved CIR-SRI map to {output_path}")

    def plot_regional_map(self, regional_data: Dict, output_path: str):
        """
        Create a simplified regional map showing WHO regions.

        Args:
            regional_data: Dictionary of regional statistics
            output_path: Output file path
        """
        fig, ax = plt.subplots(1, 1, figsize=(14, 8))

        # Simplified regional visualization
        regions = ['AFRO', 'AMRO', 'SEARO', 'WPRO', 'EURO', 'EMRO']
        colors = self.style.get_color_palette()

        # Create approximate region boxes (simplified world view)
        region_bounds = {
            'AFRO': (-20, -35, 60, 35),   # Africa
            'AMRO': (-170, -55, -30, 60), # Americas
            'SEARO': (60, -10, 145, 30),  # SE Asia
            'WPRO': (100, -45, 180, 55),  # Western Pacific
            'EURO': (-15, 35, 60, 70),    # Europe
            'EMRO': (-20, 10, 75, 40),    # Eastern Mediterranean
        }

        region_stats = {r['region']: r for r in regional_data.values()}

        for region in regions:
            bounds = region_bounds.get(region, (0, 0, 0, 0))
            stats = region_stats.get(region, {})

            mean_score = stats.get('mean_cir_sri', 0) if isinstance(stats, dict) else 0

            # Color intensity based on score
            color = colors.get(region, '#808080')

            rect = mpatches.Rectangle(
                (bounds[0], bounds[1]),
                bounds[2] - bounds[0],
                bounds[3] - bounds[1],
                linewidth=2,
                edgecolor=color,
                facecolor=color,
                alpha=0.3 + 0.5 * (mean_score / 10),
                label=f"{region}: {mean_score:.1f}"
            )
            ax.add_patch(rect)

            # Add region label
            center_x = (bounds[0] + bounds[2]) / 2
            center_y = (bounds[1] + bounds[3]) / 2
            ax.text(center_x, center_y, f"{region}\n{mean_score:.1f}",
                   ha='center', va='center', fontsize=10, fontweight='bold',
                   color='black', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        ax.set_xlim(-180, 180)
        ax.set_ylim(-60, 80)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_title('WHO Regions: Mean CIR-SRI Scores', fontsize=self.style.title_size + 2, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Legend
        legend_elements = [
            Line2D([0], [0], color=colors[region], linewidth=3, label=region)
            for region in regions
        ]
        ax.legend(handles=legend_elements, loc='lower right', title='Region')

        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved regional map to {output_path}")


class BarChartPlotter:
    """Creates bar chart visualizations."""

    def __init__(self, style: FigureStyle):
        self.style = style
        self.colors = style.get_color_palette()

    def plot_regional_comparison(self, regional_data: Dict, output_path: str):
        """
        Create bar chart comparing regions.

        Args:
            regional_data: Dictionary of regional statistics
            output_path: Output file path
        """
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        regions = list(regional_data.keys())
        mean_scores = [regional_data[r]['mean_cir_sri'] for r in regions]
        std_scores = [regional_data[r]['std_cir_sri'] for r in regions]
        region_colors = [self.colors.get(r, '#808080') for r in regions]

        # Left: Mean CIR-SRI by region
        ax = axes[0]
        bars = ax.bar(regions, mean_scores, yerr=std_scores, capsize=5,
                     color=region_colors, edgecolor='black', linewidth=0.5)

        ax.set_ylabel('Mean CIR-SRI Score')
        ax.set_xlabel('WHO Region')
        ax.set_title('Mean CIR-SRI Score by Region', fontweight='bold')
        ax.set_ylim(0, 10)
        ax.axhline(y=5, color='gray', linestyle='--', alpha=0.5, label='Midpoint')
        ax.grid(axis='y', alpha=0.3)

        # Add value labels on bars
        for bar, score in zip(bars, mean_scores):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                   f'{score:.1f}', ha='center', va='bottom', fontsize=9)

        # Right: Country count by region
        ax = axes[1]
        n_countries = [regional_data[r]['n_countries'] for r in regions]
        total_countries = [regional_data[r]['n_total_countries'] for r in regions]

        x = np.arange(len(regions))
        width = 0.35

        bars1 = ax.bar(x - width/2, n_countries, width, label='With Data',
                      color=region_colors, alpha=0.8)
        bars2 = ax.bar(x + width/2, total_countries, width, label='Total',
                      color=region_colors, alpha=0.4, hatch='//')

        ax.set_ylabel('Number of Countries')
        ax.set_xlabel('WHO Region')
        ax.set_title('Country Coverage by Region', fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(regions)
        ax.legend()
        ax.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved regional comparison to {output_path}")

    def plot_priority_tiers(self, data: pd.DataFrame, output_path: str):
        """
        Create bar chart of priority tier distribution.

        Args:
            data: DataFrame with 'region' and 'priority_tier' columns
            output_path: Output file path
        """
        fig, ax = plt.subplots(figsize=(12, 6))

        # Count by region and tier
        tier_order = ['urgent', 'high', 'medium', 'low']
        tier_colors = {'urgent': '#d73027', 'high': '#fc8d59',
                      'medium': '#fee08b', 'low': '#1a9850'}
        tier_labels = {'urgent': 'Urgent (0-3)', 'high': 'High (3-5)',
                      'medium': 'Medium (5-7)', 'low': 'Low (7-10)'}

        regions = sorted(data['region'].unique())
        x = np.arange(len(regions))
        width = 0.2

        for i, tier in enumerate(tier_order):
            counts = []
            for region in regions:
                count = len(data[(data['region'] == region) & (data['priority_tier'] == tier)])
                counts.append(count)

            bars = ax.bar(x + (i - 1.5) * width, counts, width,
                         label=tier_labels[tier], color=tier_colors[tier])

        ax.set_ylabel('Number of Countries')
        ax.set_xlabel('WHO Region')
        ax.set_title('Priority Tier Distribution by Region', fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(regions)
        ax.legend(title='Priority Tier', loc='upper right')
        ax.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved priority tier chart to {output_path}")

    def plot_dimension_comparison(self, data: pd.DataFrame, output_path: str):
        """
        Create grouped bar chart comparing dimensions.

        Args:
            data: DataFrame with dimension scores
            output_path: Output file path
        """
        fig, ax = plt.subplots(figsize=(14, 6))

        dimensions = ['sequence_availability', 'metadata_completeness',
                     'ram_site_observability', 'temporal_coverage', 'subtype_diversity']
        dim_labels = ['Seq. Avail.', 'Metadata', 'RAM Sites', 'Temporal', 'Subtypes']

        regions = sorted(data['region'].unique())
        x = np.arange(len(regions))
        width = 0.15

        dim_colors = plt.cm.Set2(np.linspace(0, 1, len(dimensions)))

        for i, (dim, label) in enumerate(zip(dimensions, dim_labels)):
            if dim in data.columns:
                scores = [data[data['region'] == r][dim].mean() for r in regions]
                ax.bar(x + (i - 2) * width, scores, width, label=label, color=dim_colors[i])

        ax.set_ylabel('Mean Score (0-2)')
        ax.set_xlabel('WHO Region')
        ax.set_title('CIR-SRI Dimension Scores by Region', fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(regions)
        ax.legend(title='Dimension', loc='upper right', ncol=3)
        ax.set_ylim(0, 2.5)
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
        ax.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved dimension comparison to {output_path}")


class HeatmapPlotter:
    """Creates heatmap visualizations."""

    def __init__(self, style: FigureStyle):
        self.style = style

    def plot_site_coverage_heatmap(self, data: pd.DataFrame, output_path: str,
                                   title: str = "RAM Site Observability by Country"):
        """
        Create heatmap of RAM site coverage.

        Args:
            data: DataFrame with countries as rows, sites as columns
            output_path: Output file path
            title: Plot title
        """
        fig, ax = plt.subplots(figsize=(12, max(8, len(data) * 0.3)))

        # Sort by total coverage
        if 'total' not in data.columns:
            data['total'] = data.sum(axis=1)
        data = data.sort_values('total', ascending=False).drop('total', axis=1)

        # Create heatmap
        cmap = LinearSegmentedColormap.from_list(
            'coverage_cmap', ['#d73027', '#fee08b', '#1a9850']
        )

        im = ax.imshow(data.values, cmap=cmap, aspect='auto', vmin=0, vmax=100)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Observability (%)', rotation=270, labelpad=15)

        # Set ticks
        ax.set_xticks(np.arange(len(data.columns)))
        ax.set_yticks(np.arange(len(data.index)))
        ax.set_xticklabels(data.columns, rotation=45, ha='right')
        ax.set_yticklabels(data.index)

        # Add value annotations
        for i in range(len(data.index)):
            for j in range(len(data.columns)):
                value = data.iloc[i, j]
                color = 'white' if value < 50 else 'black'
                ax.text(j, i, f'{value:.0f}', ha='center', va='center',
                       color=color, fontsize=7)

        ax.set_title(title, fontweight='bold', pad=10)
        ax.set_xlabel('RAM Site')
        ax.set_ylabel('Country')

        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved site coverage heatmap to {output_path}")

    def plot_temporal_heatmap(self, data: pd.DataFrame, output_path: str,
                              title: str = "Sequence Availability Over Time"):
        """
        Create heatmap of sequences by year.

        Args:
            data: DataFrame with years as columns
            output_path: Output file path
            title: Plot title
        """
        fig, ax = plt.subplots(figsize=(16, max(8, len(data) * 0.25)))

        # Sort by most recent coverage
        if len(data.columns) > 0:
            recent_years = sorted([c for c in data.columns if str(c).isdigit()])[-10:]
            data_subset = data[[c for c in recent_years if c in data.columns]]
            if 'total' in data_subset.columns:
                data_subset = data_subset.drop('total', axis=1)
            data_subset['total'] = data_subset.sum(axis=1)
            data_subset = data_subset.sort_values('total', ascending=False).drop('total', axis=1)
        else:
            return

        # Normalize for visualization
        data_norm = data_subset.apply(lambda x: np.log10(x + 1), axis=1)

        cmap = LinearSegmentedColormap.from_list(
            'temporal_cmap', ['#f7f7f7', '#fc8d59', '#d73027']
        )

        im = ax.imshow(data_norm.values, cmap=cmap, aspect='auto')

        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('log10(sequences + 1)', rotation=270, labelpad=15)

        ax.set_xticks(np.arange(len(data_norm.columns)))
        ax.set_yticks(np.arange(len(data_norm.index)))
        ax.set_xticklabels(data_norm.columns)
        ax.set_yticklabels(data_norm.index, fontsize=7)

        ax.set_title(title, fontweight='bold', pad=10)
        ax.set_xlabel('Year')
        ax.set_ylabel('Country')

        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved temporal heatmap to {output_path}")


class RadarChartPlotter:
    """Creates radar chart visualizations for regional profiles."""

    def __init__(self, style: FigureStyle):
        self.style = style
        self.colors = style.get_color_palette()

    def plot_regional_profiles(self, regional_data: Dict, output_path: str):
        """
        Create radar chart comparing regional profiles.

        Args:
            regional_data: Dictionary of regional statistics
            output_path: Output file path
        """
        dimensions = ['sequence_availability', 'metadata_completeness',
                     'ram_site_observability', 'temporal_coverage', 'subtype_diversity']
        dim_labels = ['Seq. Avail.', 'Metadata', 'RAM Sites', 'Temporal', 'Subtypes']
        n_dims = len(dimensions)

        # Create angles for radar chart
        angles = np.linspace(0, 2 * np.pi, n_dims, endpoint=False).tolist()
        angles += angles[:1]  # Complete the circle

        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))

        regions = list(regional_data.keys())

        for region in regions:
            dim_avgs = regional_data[region].get('dimension_averages', {})
            values = [dim_avgs.get(d, 0) for d in dimensions]
            values += values[:1]  # Complete the circle

            color = self.colors.get(region, '#808080')
            ax.plot(angles, values, 'o-', linewidth=2, label=region, color=color)
            ax.fill(angles, values, alpha=0.15, color=color)

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(dim_labels, size=10)
        ax.set_ylim(0, 2)
        ax.set_yticks([0.5, 1.0, 1.5, 2.0])
        ax.set_yticklabels(['0.5', '1.0', '1.5', '2.0'], size=8)
        ax.grid(color='gray', alpha=0.3)

        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0), title='Region')

        plt.title('Regional CIR-SRI Dimension Profiles', fontweight='bold', pad=20)
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved regional radar chart to {output_path}")

    def plot_country_profile(self, country_data: Dict, output_path: str):
        """
        Create radar chart for a single country.

        Args:
            country_data: Dictionary with dimension scores
            output_path: Output file path
        """
        dimensions = ['sequence_availability', 'metadata_completeness',
                     'ram_site_observability', 'temporal_coverage', 'subtype_diversity']
        dim_labels = ['Seq. Avail.', 'Metadata', 'RAM Sites', 'Temporal', 'Subtypes']
        n_dims = len(dimensions)

        angles = np.linspace(0, 2 * np.pi, n_dims, endpoint=False).tolist()
        angles += angles[:1]

        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

        values = [country_data.get(d, 0) for d in dimensions]
        values += values[:1]

        ax.plot(angles, values, 'o-', linewidth=2, color='#1976D2')
        ax.fill(angles, values, alpha=0.3, color='#1976D2')

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(dim_labels, size=10)
        ax.set_ylim(0, 2)
        ax.grid(color='gray', alpha=0.3)

        country = country_data.get('country', 'Unknown')
        plt.title(f'{country}: CIR-SRI Profile', fontweight='bold', pad=20)
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved country radar chart to {output_path}")


class DistributionPlotter:
    """Creates distribution visualizations."""

    def __init__(self, style: FigureStyle):
        self.style = style

    def plot_score_distribution(self, data: pd.DataFrame, output_path: str):
        """
        Create histogram of CIR-SRI scores.

        Args:
            data: DataFrame with 'cir_sri_score' column
            output_path: Output file path
        """
        fig, ax = plt.subplots(figsize=(10, 6))

        scores = data['cir_sri_score'].dropna()

        # Create histogram
        bins = np.arange(0, 11, 0.5)
        n, bins_edges, patches = ax.hist(scores, bins=bins, edgecolor='black',
                                          linewidth=0.5, alpha=0.7)

        # Color bars by tier
        tier_colors = {'urgent': '#d73027', 'high': '#fc8d59',
                      'medium': '#fee08b', 'low': '#1a9850'}

        for i, patch in enumerate(patches):
            bin_center = (bins_edges[i] + bins_edges[i+1]) / 2
            for tier, (min_s, max_s) in [
                ('urgent', (0, 3)), ('high', (3, 5)),
                ('medium', (5, 7)), ('low', (7, 11))
            ]:
                if min_s <= bin_center < max_s:
                    patch.set_facecolor(tier_colors[tier])
                    break

        ax.axvline(x=scores.mean(), color='red', linestyle='--',
                  linewidth=2, label=f'Mean: {scores.mean():.2f}')
        ax.axvline(x=scores.median(), color='blue', linestyle=':',
                  linewidth=2, label=f'Median: {scores.median():.2f}')

        ax.set_xlabel('CIR-SRI Score')
        ax.set_ylabel('Number of Countries')
        ax.set_title('Distribution of CIR-SRI Scores', fontweight='bold')
        ax.set_xlim(0, 10)
        ax.legend()
        ax.grid(axis='y', alpha=0.3)

        # Add tier labels
        tier_labels = [
            ('Urgent (0-3)', 1.5, '#d73027'),
            ('High (3-5)', 4, '#fc8d59'),
            ('Medium (5-7)', 6, '#fee08b'),
            ('Low (7-10)', 8.5, '#1a9850'),
        ]
        for label, x, color in tier_labels:
            ax.text(x, ax.get_ylim()[1] * 0.95, label, ha='center',
                   fontsize=8, color=color, fontweight='bold')

        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved score distribution to {output_path}")

    def plot_boxplot_by_region(self, data: pd.DataFrame, output_path: str):
        """
        Create boxplot of CIR-SRI scores by region.

        Args:
            data: DataFrame with 'cir_sri_score' and 'region' columns
            output_path: Output file path
        """
        fig, ax = plt.subplots(figsize=(12, 6))

        regions = sorted(data['region'].unique())
        colors = self.style.get_color_palette()

        box_data = [data[data['region'] == r]['cir_sri_score'].dropna().values
                   for r in regions]
        box_colors = [colors.get(r, '#808080') for r in regions]

        bp = ax.boxplot(box_data, labels=regions, patch_artist=True)

        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        ax.set_ylabel('CIR-SRI Score')
        ax.set_xlabel('WHO Region')
        ax.set_title('CIR-SRI Score Distribution by Region', fontweight='bold')
        ax.set_ylim(0, 10)
        ax.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=self.style.dpi)
        plt.close()

        logger.info(f"Saved boxplot to {output_path}")


def generate_all_figures(cir_sri_df: pd.DataFrame,
                       regional_summary: Dict,
                       site_heatmap: pd.DataFrame,
                       temporal_data: pd.DataFrame,
                       output_dir: str):
    """
    Generate all figures for the analysis.

    Args:
        cir_sri_df: CIR-SRI results DataFrame
        regional_summary: Regional summary dictionary
        site_heatmap: RAM site coverage heatmap data
        temporal_data: Temporal coverage data
        output_dir: Output directory
    """
    style = FigureStyle()
    os.makedirs(output_dir, exist_ok=True)

    # Initialize plotters
    map_plotter = WorldMapPlotter(style)
    bar_plotter = BarChartPlotter(style)
    heatmap_plotter = HeatmapPlotter(style)
    radar_plotter = RadarChartPlotter(style)
    dist_plotter = DistributionPlotter(style)

    logger.info("Generating all figures...")

    # World map
    try:
        map_plotter.plot_cir_sri_map(
            cir_sri_df, os.path.join(output_dir, 'fig1_cir_sri_world_map.png'),
            title="Global CIR-SRI Scores for HIV Capsid Inhibitor Surveillance"
        )
    except Exception as e:
        logger.warning(f"Could not generate world map: {e}")

    map_plotter.plot_regional_map(
        regional_summary, os.path.join(output_dir, 'figS1_regional_map.png')
    )

    # Bar charts
    bar_plotter.plot_regional_comparison(
        regional_summary, os.path.join(output_dir, 'fig2_regional_comparison.png')
    )

    bar_plotter.plot_priority_tiers(
        cir_sri_df, os.path.join(output_dir, 'fig3_priority_tiers.png')
    )

    bar_plotter.plot_dimension_comparison(
        cir_sri_df, os.path.join(output_dir, 'figS2_dimension_comparison.png')
    )

    # Heatmaps
    if not site_heatmap.empty:
        heatmap_plotter.plot_site_coverage_heatmap(
            site_heatmap, os.path.join(output_dir, 'fig4_ram_site_heatmap.png')
        )

    if not temporal_data.empty:
        heatmap_plotter.plot_temporal_heatmap(
            temporal_data, os.path.join(output_dir, 'figS3_temporal_coverage.png')
        )

    # Radar charts
    radar_plotter.plot_regional_profiles(
        regional_summary, os.path.join(output_dir, 'fig5_regional_profiles.png')
    )

    # Distribution plots
    dist_plotter.plot_score_distribution(
        cir_sri_df, os.path.join(output_dir, 'figS4_score_distribution.png')
    )

    dist_plotter.plot_boxplot_by_region(
        cir_sri_df, os.path.join(output_dir, 'fig6_boxplot_by_region.png')
    )

    logger.info(f"All figures saved to {output_dir}")


def main():
    """Example usage."""
    import argparse

    parser = argparse.ArgumentParser(description='Generate figures')
    parser.add_argument('cir_sri_csv', help='CIR-SRI scores CSV')
    parser.add_argument('-o', '--output', default='figures', help='Output directory')

    args = parser.parse_args()

    # Load data
    cir_sri_df = pd.read_csv(args.cir_sri_csv)

    # Load regional summary
    regional_df = pd.read_csv(args.output + '/regional_summary.csv')
    regional_summary = regional_df.set_index('region').to_dict('index')

    # Load heatmap data
    heatmap_path = args.output + '/country_heatmap.csv'
    site_heatmap = pd.read_csv(heatmap_path, index_col=0) if os.path.exists(heatmap_path) else pd.DataFrame()

    # Generate figures
    generate_all_figures(
        cir_sri_df=cir_sri_df,
        regional_summary=regional_summary,
        site_heatmap=site_heatmap,
        temporal_data=pd.DataFrame(),
        output_dir=args.output
    )

    print(f"Figures generated in {args.output}")


if __name__ == '__main__':
    main()