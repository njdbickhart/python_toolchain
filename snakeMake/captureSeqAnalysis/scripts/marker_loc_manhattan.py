"""
Modified further from Job's edits so that it works for my reporting script
"""

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

import plotly.graph_objects as go

SUGGESTIVE_LINE_LABEL = "suggestive line"
GENOMEWIDE_LINE_LABEL = "genomewide line"


def _get_hover_text(df, chrname=None, posname=None, snpname=None, genename=None, annotationname=None):
    """Format the hover text used in Manhattan and Volcano plots.
    :param (dataFrame) df: A pandas dataframe.
    :param (string) snpname: A string denoting the column name for the SNP
    names (e.g., rs number). More generally, this column could be anything
    that identifies each point being plotted. For example,
    in an Epigenomewide association study (EWAS), this could be the probe
    name or cg number. This column should be a character. This argument is
    optional, however it is necessary to specify it if you want to
    highlight points on the plot using the highlight argument in the
    figure method.
    :param (string) genename: A string denoting the column name for the
    GENE names.
    :param (string) annotationname: A string denoting the column name for
    annotations. This could be any annotation information that you
    want to include in the plot (e.g., zscore, effect size, minor allele
    frequency).
    """
    hover_text = ''
    if chrname is not None and chrname in df.columns:
        hover_text = hover_text + 'POS:' + df[chrname].astype(str) + ":" + df[posname].astype(str)

    if snpname is not None and snpname in df.columns:
        hover_text = hover_text + '<br>Sample: ' + df[snpname].astype(str)

    if genename is not None and genename in df.columns:
        hover_text = hover_text \
                     + '<br>GENE: ' \
                     + df[genename].astype(str)

    if annotationname is not None and annotationname in df.columns:
        hover_text = hover_text \
                     + '<br>' \
                     + df[annotationname].astype(str)

    return hover_text


"""
    mh = _ManhattanPlot(
        dataframe,
        chrm=chrm,
        bp=bp,
        p=p,
        snp=snp,
        gene=gene,
        annotation=annotation,
        logp=logp
    )

    return mh.figure(
        title=title,
        showgrid=showgrid,
        xlabel=xlabel,
        ylabel=ylabel,
        point_size=point_size,
        showlegend=showlegend,
        col=col,
        suggestiveline_value=suggestiveline_value,
        suggestiveline_color=suggestiveline_color,
        suggestiveline_width=suggestiveline_width,
        genomewideline_value=genomewideline_value,
        genomewideline_color=genomewideline_color,
        genomewideline_width=genomewideline_width,
        highlight=highlight,
        highlight_color=highlight_color
    )
"""


class ManhattanPlot:

    def __init__(
            self,
            dataframe,
            chrm="CHR",
            bp="BP",
            locs="LOCS",
            snp="SNP",
            gene="GENE",
            annotation=None,
            logp=True
    ):
        """
        Keyword arguments:
        - dataframe (dataframe; required): A pandas dataframe which
        must contain at least the following three columns:
            - the chromosome number
            - genomic base-pair position
            - a numeric quantity to plot such as a p-value or zscore
        - chrm (string; default 'CHR'): A string denoting the column name for the
        chromosome.  This column must be float or integer.  Minimum number
        of chromosomes required is 1. If you have X, Y, or MT chromosomes,
        be sure to renumber these 23, 24, 25, etc.
        - bp (string; default 'BP'): A string denoting the column name for the
        chromosomal position.
        - p (string; default 'P'): A string denoting the column name for the
        float quantity to be plotted on the y-axis. This column must be
        numeric. This does not have to be a p-value. It can be any
        numeric quantity such as peak heights, bayes factors, test
        statistics. If it is not a p-value, make sure to set logp = FALSE.
        - snp (string; default 'SNP'): A string denoting the column name for the
        SNP names (e.g. rs number). More generally, this column could be
        anything that identifies each point being plotted. For example, in
        an Epigenomewide association study (EWAS) this could be the probe
        name or cg number. This column should be a character. This
        argument is optional, however it is necessary to specify if you
        want to highlight points on the plot using the highlight argument
        in the figure method.
        - gene (string; default 'GENE'): A string denoting the column name for the
        GENE names. This column could be a string or a float. More
        generally, it could be any annotation information that you want
        to include in the plot.
        - annotation (string; optional): A string denoting the column name for
        an annotation. This column could be a string or a float.  This
        could be any annotation information that you want to include in
        the plot (e.g. zscore, effect size, minor allele frequency).
        - logp (bool; default True): If True, the -log10 of the p-value is
        plotted.  It isn't very useful to plot raw p-values; however,
        plotting the raw value could be useful for other genome-wide plots
        (e.g., peak heights, Bayes factors, test statistics, other
        "scores", etc.).

        Returns:
        - A ManhattanPlot object."""

        # checking the validity of the arguments

        # Make sure you have chrm, bp and p columns and that they are of
        # numeric type
        if chrm not in dataframe.columns.values:
            raise KeyError("Column %s not found in 'x' data.frame" % chrm)
        #else:
        #    if not is_numeric_dtype(dataframe[chrm].dtype):
        #        raise TypeError("%s column should be numeric. Do you have "
        #                        "'X', 'Y', 'MT', etc? If so change to "
        #                        "numbers and try again." % chrm)

        if bp not in dataframe.columns.values:
            raise KeyError("Column %s not found in 'x' data.frame" % bp)
        else:
            if not is_numeric_dtype(dataframe[bp].dtype):
                raise TypeError("%s column should be numeric type" % bp)

        if locs not in dataframe.columns.values:
            raise KeyError("Column %s not found in 'x' data.frame" % locs)
        else:
            if not is_numeric_dtype(dataframe[locs].dtype):
                raise TypeError("%s column should be numeric type" % locs)

        # Create a new DataFrame with columns named after chrm, bp, and p.
        self.data = pd.DataFrame(data=dataframe[[chrm, bp, locs]])

        if snp is not None:
            if snp not in dataframe.columns.values:
                # Warn if you don't have a snp column
                raise KeyError(
                    "snp argument specified as %s but column not found in "
                    "'x' data.frame" % snp)
            else:
                # If the input DataFrame has a snp column, add it to the new
                # DataFrame
                self.data[snp] = dataframe[snp]

        if gene is not None:
            if gene not in dataframe.columns.values:
                # Warn if you don't have a gene column
                raise KeyError(
                    "gene argument specified as %s but column not found in "
                    "'x' data.frame" % gene)
            else:
                # If the input DataFrame has a gene column, add it to the new
                # DataFrame
                self.data[gene] = dataframe[gene]

        if annotation is not None:
            if annotation not in dataframe.columns.values:
                # Warn if you don't have an annotation column
                raise KeyError(
                    "annotation argument specified as %s but column not "
                    "found in 'x' data.frame" % annotation
                )
            else:
                # If the input DataFrame has a gene column, add it to the new
                # DataFrame
                self.data[annotation] = dataframe[annotation]

        self.xlabel = ""
        self.ticks = []
        self.ticksLabels = []
        self.nChr = len(dataframe[chrm].unique())
        self.chrName = chrm
        self.lName = locs
        self.snpName = snp
        self.geneName = gene
        self.annotationName = annotation
        self.logp = logp

        # Set positions, ticks, and labels for plotting

        self.index = 'INDEX'
        self.pos = 'POSITION'

        # Fixes the bug where one chromosome is missing by adding a sequential
        # index column.
        idx = 0
        for i in self.data[chrm].unique():
            idx = idx + 1
            self.data.loc[self.data[chrm] == i, self.index] = int(idx)
        # Set the type to be the same as provided for chrm column
        self.data[self.index] = \
            self.data[self.index].astype(int)

        # This section sets up positions and ticks. Ticks should be placed in
        # the middle of a chromosome. The new pos column is added that keeps
        # a running sum of the positions of each successive chromosome.
        # For example:
        # chrm bp pos
        # 1   1  1
        # 1   2  2
        # 2   1  3
        # 2   2  4
        # 3   1  5

        if self.nChr == 1:
            # For a single chromosome
            self.data[self.pos] = self.data[bp]
            self.ticks.append(int(len(self.data[self.pos]) / 2.) + 1)
            self.xlabel = "Chromosome %s position" % (self.data[chrm].unique())
            self.ticksLabels = self.ticks
        else:
            # For multiple chromosomes
            lastbase = 0
            for i in self.data[self.index].unique():
                if i == 1:
                    self.data.loc[self.data[self.index] == i, self.pos] = \
                        self.data.loc[self.data[self.index] == i, bp].values
                else:
                    prevbp = self.data.loc[self.data[self.index] == i - 1, bp]
                    # Shift the basepair position by the largest bp of the
                    # current chromosome
                    lastbase = lastbase + prevbp.iat[-1]

                    self.data.loc[self.data[self.index] == i, self.pos] = \
                        self.data.loc[self.data[self.index] == i, bp].values \
                        + lastbase

                tmin = min(self.data.loc[self.data[self.index] == i, self.pos])
                tmax = max(self.data.loc[self.data[self.index] == i, self.pos])
                self.ticks.append(int((tmin + tmax) / 2.) + 1)

            self.xlabel = 'Chromosome'
            self.data[self.pos] = self.data[self.pos].astype(
                self.data[bp].dtype)

            chromList = list(self.data[chrm].unique())
            #if self.nChr > 10:  # To avoid crowded labels
            #    self.ticksLabels = [
            #        t if np.mod(int(i+ 1), 2) or int(i +1) < 6  # Only every two ticks
            #        else ''
            #        for i, t in enumerate(self.data[chrm].unique())
            #    ]
            #else:
            self.ticksLabels = self.data[chrm].unique()  # All the ticks

    def update_p_col(self, l_col):
        """Keyword arguments:
        - alternate_p_col (string, default ""): different column than the
            standard to use to generate a plot.

        :return: None, updates internal p_col
        """
        self.lName = l_col

    def get_traces(
            self,
            point_size=5,
            showlegend=True,
            col=None,
            genomewideline_value=-np.log(1),
            highlight=True,
            highlight_color="red"
    ):
        """Keyword arguments:
    - point_size (number; default 5): Size of the points of the
        scatter plot.
    - showlegend (bool; default True): Boolean indicating whether
        legends should be shown.
    - col (string; optional): A string representing the color of the
        points of the Scatter plot. Can be in any color format
        accepted by plotly.graph_objects.
    - genomewideline_value (bool | float; default -log10(5e-8)): A
        boolean which must be either False to deactivate the option, or a
        numerical value corresponding to the p-value above which the
        data points are considered significant.
    - highlight (bool; default True): Whether to turn on or off the
        highlighting of data points considered significant.
    - highlight_color (string; default 'red'): Color of the data
        points highlighted because they are significant. Can be in any
        color format accepted by plotly.graph_objects.

    Returns:
    - A list of traces to be used to create a Manhattan plot.

        """

        data_to_plot = []  # To contain the data traces
        tmp = pd.DataFrame()  # Empty DataFrame to contain the highlighted data

        if highlight:
            if not isinstance(highlight, bool):
                if self.snpName not in self.data.columns.values:
                    raise KeyError(
                        "snp argument specified for highlight as %s but "
                        "column not found in the data.frame" % self.snpName
                    )
            else:
                #if not genomewideline_value:
                #    raise Warning(
                #        "The genomewideline_value you entered is not a "
                #        "positive value, or False, you cannot set highlight "
                #        "to True in that case.")
                tmp = self.data

                # Sort the p-values (or -log10(p-values) above the line
                if genomewideline_value:
                    if self.logp:
                        tmp = tmp.loc[-np.log(tmp[self.lName])
                                      >= genomewideline_value]
                    else:
                        tmp = tmp.loc[tmp[self.lName] >= genomewideline_value]

                highlight_hover_text = _get_hover_text(
                    tmp,
                    chrname=self.chrName, 
                    posname=self.pos,
                    snpname=self.snpName,
                    genename=self.geneName,
                    annotationname=self.annotationName
                )

                if not tmp.empty:
                    data_to_plot.append(
                        go.Scattergl(
                            x=tmp[self.pos].values,
                            y=-np.log(tmp[self.lName].values) if self.logp
                            else tmp[self.lName].values,
                            mode="markers",
                            showlegend=showlegend,
                            text=highlight_hover_text,
                            marker=dict(
                                color=highlight_color,
                                size=point_size
                            ),
                            name="Point(s) of interest"
                        )
                    )

        # Remove the highlighted data from the DataFrame if not empty
        if tmp.empty:
            data = self.data
        else:
            if self.logp:
                data = self.data.loc[-np.log(self.data[self.lName])
                                     < genomewideline_value]
            else:
                data = self.data.loc[self.data[self.lName]
                                     < genomewideline_value]

        if self.nChr == 1:

            if col is None:
                col = ['black']

            hover_text = _get_hover_text(
                data,
                chrname=self.chrName, 
                posname=self.pos,
                snpname=self.snpName,
                genename=self.geneName,
                annotationname=self.annotationName
            )

            data_to_plot.append(
                go.Scattergl(
                    x=data[self.pos].values,
                    y=-np.log(data[self.lName].values) if self.logp
                    else data[self.lName].values,
                    mode="markers",
                    showlegend=showlegend,
                    name="chr%s" % data[self.chrName].unique(),
                    marker={
                        'color': col[0],
                        'size': point_size
                    },
                    text=hover_text
                )
            )
        else:
            icol = 0
            if col is None:
                col = [
                    'black' if np.mod(i, 2)
                    else 'grey' for i in range(self.nChr)
                ]

            for i in data[self.index].unique():
                tmp = data[data[self.index] == i]

                chromo = tmp[self.chrName].unique()  # Get chromosome name

                hover_text = _get_hover_text(
                    tmp,
                    chrname=self.chrName, 
                    posname=self.pos,
                    snpname=self.snpName,
                    genename=self.geneName,
                    annotationname=self.annotationName
                )

                data_to_plot.append(
                    go.Scattergl(
                        x=tmp[self.pos].values,
                        y=-np.log(tmp[self.lName].values) if self.logp
                        else tmp[self.lName].values,
                        mode="markers",
                        showlegend=showlegend,
                        name="Chr%s" % chromo,
                        marker={
                            'color': col[icol],
                            'size': point_size
                        },
                        text=hover_text
                    )
                )

                icol = icol + 1
        return data_to_plot

    def get_threshold_lines(
            self,
            suggestiveline_value=-np.log(2),
            suggestiveline_color='blue',
            suggestiveline_width=1,
            genomewideline_value=-np.log(1),
            genomewideline_color='red',
            genomewideline_width=1,
            axis=1,
    ) -> list:
        """Keyword arguments:
     - suggestiveline_value (bool | float; default 8): A value which
         must be either False to deactivate the option, or a numerical value
         corresponding to the p-value at which the line should be
         drawn. The line has no influence on the data points.
     - suggestiveline_color (string; default 'grey'): Color of the
         suggestive line.
     - suggestiveline_width (number; default 2): Width of the
         suggestive line.
     - genomewideline_value (bool | float; default -log10(5e-8)): A
         boolean which must be either False to deactivate the option, or a
         numerical value corresponding to the p-value above which the
         data points are considered significant.
     - genomewideline_color (string; default 'red'): Color of the
         genome-wide line. Can be in any color format accepted by
         plotly.graph_objects.
     - genomewideline_width (number; default 1): Width of the genome
       wide line.
     - axis: (number; default 1): The axis to generate layout for

     Returns:
     - A list with threshold lines.

         """

        xmin = min(self.data[self.pos].values)
        xmax = max(self.data[self.pos].values)
        xref = "x" if axis == 1 else f"x{axis}"
        yref = "y" if axis == 1 else f"y{axis}"

        horizontallines = []

        if suggestiveline_value:
            suggestiveline = go.layout.Shape(
                name=SUGGESTIVE_LINE_LABEL,
                type="line",
                fillcolor=suggestiveline_color,
                line=dict(
                    color=suggestiveline_color,
                    width=suggestiveline_width
                ),
                x0=xmin, x1=xmax, xref=xref,
                y0=suggestiveline_value, y1=suggestiveline_value, yref=yref
            )
            horizontallines.append(suggestiveline)

        if genomewideline_value:
            genomewideline = go.layout.Shape(
                name=GENOMEWIDE_LINE_LABEL,
                type="line",
                fillcolor=genomewideline_color,
                line=dict(
                    color=genomewideline_color,
                    width=genomewideline_width
                ),
                x0=xmin, x1=xmax, xref=xref,
                y0=genomewideline_value, y1=genomewideline_value, yref=yref
            )
            horizontallines.append(genomewideline)
        return horizontallines

    def get_x_layout(
            self,
            showgrid=True,
            xlabel=None,
    ) -> dict:
        """Keyword arguments:
     - title (string; default 'Manhattan Plot'): The title of the
         graph.
     - showgrid (bool; default True): Boolean indicating whether
         gridlines should be shown.
     - xlabel (string; optional): Label of the x axis.
     - ylabel (string; default '-log10(p)'): Label of the y axis.
     - axis: (number; default 1): The axis to generate layout for

     Returns:
     - A dict with layout elements, can be used to make a layout instance.

         """
        xmin = min(self.data[self.pos].values)
        xmax = max(self.data[self.pos].values)

        if self.nChr == 1:
            # If single chromosome, ticks and labels automatic.
            return {
                'title': self.xlabel if xlabel is None else xlabel,
                'showgrid': showgrid,
                'range': [xmin, xmax],
                'showticklabels': True
            }
        # if multiple chroms, use the ticks and labels created previously.
        return {
            'title': self.xlabel if xlabel is None else xlabel,
            'showgrid': showgrid,
            'range': [xmin, xmax],
            'tickmode': "array",
            'tickvals': self.ticks,
            'ticktext': self.ticksLabels,
            'ticks': "outside",
            'showticklabels': True
        }


    def get_layout(
            self,
            title="Manhattan Plot",
            showgrid=True,
            xlabel=None,
            ylabel='-log10(p)',
            suggestiveline_value=-np.log(2),
            suggestiveline_color='blue',
            suggestiveline_width=1,
            genomewideline_value=-np.log(1),
            genomewideline_color='red',
            genomewideline_width=1,
            hovermode='closest',
            axis=1,
    ) -> dict:
        """Keyword arguments:
     - title (string; default 'Manhattan Plot'): The title of the
         graph.
     - showgrid (bool; default True): Boolean indicating whether
         gridlines should be shown.
     - xlabel (string; optional): Label of the x axis.
     - ylabel (string; default '-log10(p)'): Label of the y axis.
     - suggestiveline_value (bool | float; default 8): A value which
         must be either False to deactivate the option, or a numerical value
         corresponding to the p-value at which the line should be
         drawn. The line has no influence on the data points.
     - suggestiveline_color (string; default 'grey'): Color of the
         suggestive line.
     - suggestiveline_width (number; default 2): Width of the
         suggestive line.
     - genomewideline_value (bool | float; default -log10(5e-8)): A
         boolean which must be either False to deactivate the option, or a
         numerical value corresponding to the p-value above which the
         data points are considered significant.
     - genomewideline_color (string; default 'red'): Color of the
         genome-wide line. Can be in any color format accepted by
         plotly.graph_objects.
     - genomewideline_width (number; default 1): Width of the genome
       wide line.
     - axis: (number; default 1): The axis to generate layout for

     Returns:
     - A dict with layout elements, can be used to make a layout instance.

         """
        horizontallines = self.get_threshold_lines(
            suggestiveline_value=suggestiveline_value,
            suggestiveline_color=suggestiveline_color,
            suggestiveline_width=suggestiveline_width,
            genomewideline_value=genomewideline_value,
            genomewideline_color=genomewideline_color,
            genomewideline_width=genomewideline_width,
            axis=axis
        )

        x_layout = self.get_x_layout(
            showgrid=showgrid,
            xlabel=xlabel,
        )

        xaxis = "xaxis" if axis == 1 else f"xaxis{axis}"
        yaxis = "yaxis" if axis == 1 else f"yaxis{axis}"
        return {
            "title": title,
            xaxis: x_layout,
            yaxis: {'title': ylabel},
            "hovermode": hovermode,
            "shapes": horizontallines
        }

    def figure(
            self,
            title="Manhattan Plot",
            showgrid=True,
            xlabel=None,
            ylabel='-log(Insertions)',
            point_size=5,
            showlegend=True,
            col=None,
            suggestiveline_value=-np.log(2),
            suggestiveline_color='blue',
            suggestiveline_width=1,
            genomewideline_value=-np.log(1),
            genomewideline_color='red',
            genomewideline_width=1,
            highlight=True,
            highlight_color="red",
    ):
        """Keyword arguments:
    - title (string; default 'Manhattan Plot'): The title of the
        graph.
    - showgrid (bool; default True): Boolean indicating whether
        gridlines should be shown.
    - xlabel (string; optional): Label of the x axis.
    - ylabel (string; default '-log10(p)'): Label of the y axis.
    - point_size (number; default 5): Size of the points of the
        scatter plot.
    - showlegend (bool; default True): Boolean indicating whether
        legends should be shown.
    - col (string; optional): A string representing the color of the
        points of the Scatter plot. Can be in any color format
        accepted by plotly.graph_objects.
    - suggestiveline_value (bool | float; default 8): A value which
        must be either False to deactivate the option, or a numerical value
        corresponding to the p-value at which the line should be
        drawn. The line has no influence on the data points.
    - suggestiveline_color (string; default 'grey'): Color of the
        suggestive line.
    - suggestiveline_width (number; default 2): Width of the
        suggestive line.
    - genomewideline_value (bool | float; default -log10(5e-8)): A
        boolean which must be either False to deactivate the option, or a
        numerical value corresponding to the p-value above which the
        data points are considered significant.
    - genomewideline_color (string; default 'red'): Color of the
        genome-wide line. Can be in any color format accepted by
        plotly.graph_objects.
    - genomewideline_width (number; default 1): Width of the genome
      wide line.
    - highlight (bool; default True): Whether to turn on or off the
        highlighting of data points considered significant.
    - highlight_color (string; default 'red'): Color of the data
        points highlighted because they are significant. Can be in any
        color format accepted by plotly.graph_objects.

    Returns:
    - A figure formatted for plotly.graph_objects.

        """
        data_to_plot = self.get_traces(
            point_size=point_size,
            showlegend=showlegend,
            col=col,
            genomewideline_value=genomewideline_value,
            highlight=highlight,
            highlight_color=highlight_color
        )

        layout = self.get_layout(
            title=title,
            showgrid=showgrid,
            xlabel=xlabel,
            ylabel=ylabel,
            suggestiveline_value=suggestiveline_value,
            suggestiveline_color=suggestiveline_color,
            suggestiveline_width=suggestiveline_width,
            genomewideline_value=genomewideline_value,
            genomewideline_color=genomewideline_color,
            genomewideline_width=genomewideline_width,
            hovermode='closest'
        )

        return go.Figure(data=data_to_plot, layout=go.Layout(layout))