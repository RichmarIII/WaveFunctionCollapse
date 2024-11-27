#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
#include <bitset>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <queue>
#include <string>
#include <chrono>
#include <random>
#include <fstream>
#include <optional>
#include <stack>

// Include stb_image for image loading
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Include stb_image_write for image writing
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace PSE::WFC
{

    class WaveFunctionCollapse
    {
    public:
        /**
         * @brief Represents an index of a pattern.
         *
         */
        using PatternIdx = std::size_t;

        /**
         * @brief Represents the size of a pattern in pixels.
         *
         */
        using PatternSize = std::size_t;

        /**
         * @brief Represents a pattern. A pattern is a block of RGBA pixels.
         *
         * The actual data of a pattern
         *
         */
        using PatternData = std::vector<std::tuple<std::uint8_t, std::uint8_t, std::uint8_t, std::uint8_t>>;

        /**
         * @brief Represents the edge data of a pattern.
         *
         * The edge data is stored in a tuple.
         * The tuple contains the RGBA values of the particular edge of the pattern.
         *
         */
        using PatternEdgeData = std::vector<std::tuple<std::uint8_t, std::uint8_t, std::uint8_t, std::uint8_t>>;

        /**
         * @brief Represents the coordinates of a pattern in the source image.
         *
         * The coordinates are in (x, y) format.
         *
         */
        using PatternCoord = std::tuple<std::size_t, std::size_t>;

        /**
         * @brief Represents the adjacency rules for a particular side of a pattern.
         *
         * Adjacency rules specify which patterns can be placed next to each other.
         * The key is a pattern index, and the value is a vector of pattern indices that can be placed next to the specific side of the key pattern.
         *
         */
        using SideAdjacencyMap = std::unordered_map<PatternIdx, std::vector<PatternIdx>>;

        /**
         * @brief Represents the amount of cells in a grid in a single dimension.
         *
         */
        using GridSize = std::size_t;

        /**
         * @brief Represents the actual pixel data of an image.
         *
         * The image data is stored in a tuple.
         * The first element is the width of the image.
         * The second element is the height of the image.
         * The third element is a vector of tuples, where each tuple represents an RGBA pixel.
         * The RGBA pixel is stored as (R, G, B, A) values.
         * The pixel data is stored in row-major order.
         *
         */
        using InputImage = std::tuple<std::size_t, std::size_t, std::vector<std::tuple<std::uint8_t, std::uint8_t, std::uint8_t, std::uint8_t>>>;

        enum class EEdge
        {
            Left,
            Right,
            Top,
            Bottom
        };

        /**
         * @brief Represents a view into a unique pattern.
         *
         * A pattern is a block of RGBA pixels.
         * This struct doesn't store the actual pixel data, only the indices of the pattern.
         * The actual pixel data is stored elsewhere (e.g. in a source image). This struct is used to refer to the pattern.
         *
         */
        struct UniquePatternView
        {
            /**
             * @brief The index of the pattern.
             *
             * This index is used to refer to the unique pattern in the source image.
             * Index is calculated by iterating over the source image from left to right, top to bottom.
             * The first pattern encountered is assigned index 0, the second pattern is assigned index 1, and so on.
             *
             */
            PatternIdx Idx{};

            /**
             * @brief The frequency of the pattern in the source image.
             *
             * The frequency is the number of times the pattern appears in the source image.
             *
             */
            std::size_t Frequency{};

            /**
             * @brief Calculates and returns the coordinates of the pattern in the source image.
             *
             * @param NumPatternsPerRow The number of patterns in each row of the source image.
             * @param NumPatternsPerCol The number of patterns in each column of the source image.
             *
             * @return The coordinates of the pattern in the source image.
             *
             * @example
             * // Example usage:
             * // Assume the source image has 10 patterns per row.
             * PatternView PV{5};
             * auto Coord = PV.CalcPatternCoord(10);
             * // Coord is now (5, 0)
             *
             * PatternView PV{15};
             * auto Coord = PV.CalcPatternCoord(10);
             * // Coord is now (5, 1)
             *
             */
            static PatternCoord ComputePatternCoord(const PatternIdx Idx, std::size_t NumPatternsPerRow, std::size_t NumPatternsPerCol)
            {
                assert(NumPatternsPerRow > 0 && "NumPatternsPerRow must be greater than 0");
                assert(NumPatternsPerCol > 0 && "NumPatternsPerCol must be greater than 0");
                assert(Idx < NumPatternsPerRow * NumPatternsPerCol && "Idx must be less than NumPatternsPerRow * NumPatternsPerCol");

                const auto PatternX = Idx % NumPatternsPerRow;
                const auto PatternY = Idx / NumPatternsPerRow;

                return {PatternX, PatternY};
            }

            /**
             * @brief Gets the pattern data from the input image.
             *
             * @param Coord The coordinates of the pattern in the source image.
             * @param PatSize The size of the pattern.
             * @param InputImage The input image data.
             *
             * @return PatternData The pattern data extracted from the input image.
             */
            static PatternData GetPatternData(PatternCoord Coord, const PatternSize PatSize, const InputImage &InputImage)
            {
                const auto Width = std::get<0>(InputImage);
                const auto &ImagePixels = std::get<2>(InputImage);

                const auto [PatternX, PatternY] = Coord;

                PatternData PData{};
                PData.reserve(PatSize * PatSize);

                const auto StartX = PatternX * (PatSize - 1);
                const auto StartY = PatternY * (PatSize - 1);

                for (std::size_t py = 0; py < PatSize; ++py)
                {
                    for (std::size_t px = 0; px < PatSize; ++px)
                    {
                        // Get the pixel data from the window
                        auto &Pixel = ImagePixels[(StartY + py) * Width + (StartX + px)];
                        PData.emplace_back(Pixel);
                    }
                }

                return PData;
            }

            static PatternEdgeData GetPatternEdgeData(const PatternData &Data, const PatternSize PatSize, const EEdge Edge)
            {
                PatternEdgeData EdgeData{};
                EdgeData.reserve(PatSize);

                switch (Edge)
                {
                case EEdge::Left:
                    for (std::size_t i = 0; i < PatSize; ++i)
                    {
                        EdgeData.push_back(Data[i * PatSize]);
                    }
                    break;
                case EEdge::Right:
                    for (std::size_t i = 0; i < PatSize; ++i)
                    {
                        EdgeData.push_back(Data[i * PatSize + PatSize - 1]);
                    }
                    break;
                case EEdge::Top:
                    EdgeData.insert(EdgeData.end(), Data.begin(), Data.begin() + PatSize);
                    break;
                case EEdge::Bottom:
                    EdgeData.insert(EdgeData.end(), Data.end() - PatSize, Data.end());
                    break;
                default:
                    break;
                }

                return EdgeData;
            }
        };

        /**
         * @brief Represents a cell in the output grid.
         *
         * A cell can be in one of two states:
         * 1. Collapsed: The cell has been collapsed to a single pattern.
         * 2. Uncollapsed: The cell has multiple possible patterns that can be placed in it.
         *
         */
        struct Cell
        {
            /// @brief The possible patterns that can be placed in this cell. Only valid if the cell is not collapsed.
            std::vector<PatternIdx> PossiblePatterns{};

            /// @brief The index of the pattern that was collapsed into this cell.
            std::optional<PatternIdx> CollapsedPatternIdx{};
        };

        /**
         * @brief Represents the adjacency rules for patterns.
         *
         * Adjacency rules specify which patterns can be placed next to each other.
         *
         */
        struct AdjacencyRules
        {
            /**
             * @brief The patterns that can be placed to the left of each pattern.
             *
             * The key is a pattern index, and the value is a vector of pattern indices that can be placed to the left of the key pattern.
             */
            SideAdjacencyMap Left{};

            /**
             * @brief The patterns that can be placed to the right of each pattern.
             *
             * The key is a pattern index, and the value is a vector of pattern indices that can be placed to the right of the key pattern.
             */
            SideAdjacencyMap Right{};

            /**
             * @brief The patterns that can be placed above each pattern.
             *
             * The key is a pattern index, and the value is a vector of pattern indices that can be placed above the key pattern.
             */
            SideAdjacencyMap Top{};

            /**
             * @brief The patterns that can be placed below each pattern.
             *
             * The key is a pattern index, and the value is a vector of pattern indices that can be placed below the key pattern.
             */
            SideAdjacencyMap Bottom{};
        };

        struct RunParameters
        {
            RunParameters(std::string_view InputImageFilename, std::string_view OutputImageFilename, PatternSize PatSize, GridSize GridOutputSize)
                : InputImageFilename(InputImageFilename), OutputImageFilename(OutputImageFilename), PatSize(PatSize), GridSize(GridOutputSize) {}

            std::string InputImageFilename{};
            std::string OutputImageFilename{};

            /**
             * @brief The size of the patterns to extract from the input image.
             *
             * The Wave Function Collapse algorithm extracts patterns from the input image.
             * This parameter specifies the size of the patterns to extract from the image.
             * A pattern is a block of RGBA pixels. The size of the pattern is specified in pixels.
             * A size of 2 means that each pattern is a 2x2 block of pixels.
             *
             */
            PatternSize PatSize{};

            /**
             * @brief The size of the output grid (in cells).
             *
             * The output grid is a square grid of cells. This parameter specifies the size of the grid.
             * The grid size is specified in cells. A size of 3 means that the grid is 3x3 cells in size.
             *
             */
            GridSize GridSize{};
        };

        /**
         * @brief Represents a grid of cells. This is the output of the Wave Function Collapse algorithm.
         *
         * The grid is a 2D vector of cells. Each cell can be in one of two states:
         * 1. Collapsed: The cell has been collapsed to a single pattern.
         * 2. Uncollapsed: The cell has multiple possible patterns that can be placed in it.
         *
         * The grid is initialized with all cells in the uncollapsed state.
         * The Wave Function Collapse algorithm collapses the cells one by one until the grid is fully collapsed.
         *
         * This grid is always indexed in [x][y] order (column-major order).
         *
         */
        using Grid = std::vector<std::vector<Cell>>;

        WaveFunctionCollapse(GridSize GridOutputSize);
        WaveFunctionCollapse(const WaveFunctionCollapse &) = delete;
        WaveFunctionCollapse(WaveFunctionCollapse &&) = delete;
        WaveFunctionCollapse &operator=(const WaveFunctionCollapse &) = delete;
        WaveFunctionCollapse &operator=(WaveFunctionCollapse &&) = delete;

        static bool Run(const RunParameters &Params);

    private:
        static InputImage LoadImage(std::string_view InputImageFilename);

        /**
         * @brief Extracts patterns from the input image. Each pattern is a block of RGBA pixels.
         *
         * @param Image The input image data.
         * @param PatSize The size of the patterns to extract.
         *
         * @return A vector of unique patterns extracted from the input image.
         *
         */
        static std::vector<UniquePatternView> ExtractPatterns(const InputImage &Image, PatternSize PatSize);

        /**
         * @brief Extracts and returns adjacency rules from the patterns.
         *
         * @param UniquePatterns The unique patterns extracted from the input image.
         * @param PatSize The size of the patterns.
         * @param InputImage The input image data.
         *
         * @return The adjacency rules extracted from the patterns.
         *
         */
        static AdjacencyRules ExtractAdjacencyRules(const std::vector<UniquePatternView> &UniquePatterns, const PatternSize PatSize, const WaveFunctionCollapse::InputImage &InputImage);

        /**
         * @brief Initializes the grid with all patterns as possible in each cell.
         *
         * @param GridSize The size of the output grid (in cells).
         * @param UniquePatterns The unique patterns extracted from the input image.
         *
         */
        static Grid InitializeGrid(const GridSize GridSize, const std::vector<UniquePatternView> &UniquePatterns);

        /**
         * @brief Collapses the grid by propagating constraints.
         *
         * @param Grid The grid to collapse.
         * @param AdjacencyRules The adjacency rules extracted from the patterns.
         * @param CellX The X-coordinate of the cell that was collapsed.
         * @param CellY The Y-coordinate of the cell that was collapsed.
         *
         * @return True if the grid was successfully propagated, false otherwise.
         *
         */
        static bool Propagate(Grid &Grid, const AdjacencyRules &AdjacencyRules, std::size_t CellX, std::size_t CellY);

        static bool PropagateWithBackTrack(Grid &Grid, const AdjacencyRules &AdjacencyRules, std::size_t CellX, std::size_t CellY);

        /**
         * @brief Chooses the next cell to collapse. Typically, the cell with the lowest entropy is chosen.
         *
         * @param Grid The grid to choose the next cell from.
         *
         * @return A tuple containing the X and Y coordinates of the next cell to collapse.
         */
        static std::tuple<std::size_t, std::size_t> ChooseNextCell(const Grid &Grid);

        /**
         * @brief Visualizes the entropy values of the grid.
         *
         * This function prints the entropy values of each cell in the grid.
         * The entropy value is the number of possible patterns that can be placed in the cell.
         *
         * @param Grid The grid to visualize.
         *
         */
        static void VisualizeGridEntropy(const Grid &Grid);

        /**
         * @brief Visualizes the adjacency rules.
         *
         * @param Rules The adjacency rules to visualize.
         */
        static void VisualizeAdjacencyRules(const AdjacencyRules &Rules);
    };

    WaveFunctionCollapse::InputImage WaveFunctionCollapse::LoadImage(std::string_view InputImageFilename)
    {
        // Load the input image
        int Width{}, Height{}, Channels{};
        unsigned char *ImageData = stbi_load(InputImageFilename.data(), &Width, &Height, &Channels, 4);

        if (!ImageData)
        {
            std::cerr << "Failed to load image: " << InputImageFilename << "\n";
            return {};
        }

        // Copy the image data into the input image
        InputImage Image = {std::size_t(Width), std::size_t(Height), {}};
        auto &ImagePixels = std::get<2>(Image);
        ImagePixels.reserve(Width * Height);

        for (int y = 0; y < Height; ++y)
        {
            for (int x = 0; x < Width; ++x)
            {
                auto *PixelData = ImageData + (x + y * Width) * 4; // RGBA
                ImagePixels.emplace_back(PixelData[0], PixelData[1], PixelData[2], PixelData[3]);
            }
        }

        // Free the image data
        stbi_image_free(ImageData);

        return Image;
    }

    std::vector<WaveFunctionCollapse::UniquePatternView> WaveFunctionCollapse::ExtractPatterns(const InputImage &Image, PatternSize PatSize)
    {
        // Load the input image
        const auto Width = std::get<0>(Image);
        const auto Height = std::get<1>(Image);
        auto &ImagePixels = std::get<2>(Image);

        // Extract unique patterns from the input image
        std::cout << "Extracting unique patterns...\n";

        // Used to store the unique patterns and their frequencies
        std::vector<UniquePatternView> UniquePatterns{};
        std::vector<PatternData> PatternDatas{};
        PatternIdx CurrentPatternIdx{0};

        // Extract unique patterns from the input image

        // We need to move the "Pattern Window" over the image and extract patterns from each window
        // The window size is specified by the pattern size
        // The window moves from left to right, top to bottom
        // The window moves by PatSize - 1 pixels horizontally and vertically
        // This ensures that we don't miss any patterns and only overlap the edges

        // Loop through the image with a sliding window
        for (std::size_t y = 0; y <= Height - PatSize; y += PatSize - 1)
        {
            for (std::size_t x = 0; x <= Width - PatSize; x += PatSize - 1)
            {
                // Extract the pattern from the window
                PatternData PData{};
                for (std::size_t py = 0; py < PatSize; ++py)
                {
                    for (std::size_t px = 0; px < PatSize; ++px)
                    {
                        // Get the pixel data from the window
                        auto &Pixel = ImagePixels[(y + py) * Width + (x + px)];
                        PData.emplace_back(Pixel);
                    }
                }

                // If the pattern is unique, add it to the list of unique patterns
                auto Iter = std::find(PatternDatas.begin(), PatternDatas.end(), PData);
                if (Iter == PatternDatas.end())
                {
                    UniquePatternView PatternView{CurrentPatternIdx, 1};
                    UniquePatterns.emplace_back(PatternView);
                    PatternDatas.emplace_back(PData);
                }
                else
                {
                    // If the pattern is not unique, increment the frequency of the existing pattern

                    // Get the index of the existing pattern
                    auto ExistingPatternIdx = std::distance(PatternDatas.begin(), Iter);

                    // Increment the frequency of the existing pattern
                    UniquePatterns[ExistingPatternIdx].Frequency++;
                }

                // Increment the current pattern index
                ++CurrentPatternIdx;
            }
        }

        std::cout << "Number of unique patterns after deduplication: " << UniquePatterns.size() << "\n";

        // Verify that there are no duplicate patterns
        for (std::size_t i = 0; i < UniquePatterns.size(); ++i)
        {
            for (std::size_t j = i + 1; j < UniquePatterns.size(); ++j)
            {
                if (PatternDatas[i] == PatternDatas[j])
                {
                    std::cerr << "Duplicate pattern found at indices " << i << " and " << j << "\n";
                }
            }
        }

        return UniquePatterns;
    }

    WaveFunctionCollapse::AdjacencyRules WaveFunctionCollapse::ExtractAdjacencyRules(const std::vector<UniquePatternView> &UniquePatterns, const PatternSize PatSize, const WaveFunctionCollapse::InputImage &InputImage)
    {
        // Extract adjacency rules from the patterns
        std::cout << "Extracting adjacency rules...\n";

        // Initialize the adjacency rules
        AdjacencyRules Rules{};
        // Rules.Left.reserve(UniquePatterns.size());  // Don't reserve, could be extremely large
        // Rules.Right.reserve(UniquePatterns.size()); // Don't reserve, could be extremely large
        // Rules.Top.reserve(UniquePatterns.size()); // Don't reserve, could be extremely large
        // Rules.Bottom.reserve(UniquePatterns.size()); // Don't reserve, could be extremely large

        // Loop through each pattern and extract adjacency rules.
        // This is not the unique patterns, but the actual patterns.
        auto &ImagePixels = std::get<2>(InputImage);
        const auto Width = std::get<0>(InputImage);
        const auto Height = std::get<1>(InputImage);
        const auto NumPatternsPerRow = ((Width - PatSize) / (PatSize - 1)) + 1;
        const auto NumPatternsPerCol = ((Height - PatSize) / (PatSize - 1)) + 1;

        for (PatternIdx i = 0; i < UniquePatterns.size(); ++i)
        {
            auto &UniquePatternOuter = UniquePatterns[i];
            auto CurrentCoord = UniquePatternView::ComputePatternCoord(UniquePatternOuter.Idx, NumPatternsPerRow, NumPatternsPerCol);

            // Extract the pattern data
            auto PatternDataOuter = UniquePatternView::GetPatternData(CurrentCoord, PatSize, InputImage);

            for (PatternIdx j = 0; j < UniquePatterns.size(); ++j)
            {
                auto &UniquePatternInner = UniquePatterns[j];
                auto PatternDataInner = UniquePatternView::GetPatternData(UniquePatternView::ComputePatternCoord(UniquePatternInner.Idx, NumPatternsPerRow, NumPatternsPerCol), PatSize, InputImage);

                if (UniquePatternInner.Idx == UniquePatternOuter.Idx)
                {
                    // Skip the same pattern
                    continue;
                }

                // Check if the inner pattern can be placed to the left of the outer pattern by comparing the right side of the inner pattern with the left side of the outer pattern
                auto InnerRightEdge = UniquePatternView::GetPatternEdgeData(PatternDataInner, PatSize, EEdge::Right);
                auto OuterLeftEdge = UniquePatternView::GetPatternEdgeData(PatternDataOuter, PatSize, EEdge::Left);
                if (InnerRightEdge == OuterLeftEdge)
                {
                    Rules.Left[UniquePatternOuter.Idx].push_back(UniquePatternInner.Idx);
                }

                // Check if the inner pattern can be placed to the right of the outer pattern by comparing the left side of the inner pattern with the right side of the outer pattern
                auto InnerLeftEdge = UniquePatternView::GetPatternEdgeData(PatternDataInner, PatSize, EEdge::Left);
                auto OuterRightEdge = UniquePatternView::GetPatternEdgeData(PatternDataOuter, PatSize, EEdge::Right);
                if (InnerLeftEdge == OuterRightEdge)
                {
                    Rules.Right[UniquePatternOuter.Idx].push_back(UniquePatternInner.Idx);
                }

                // Check if the inner pattern can be placed above the outer pattern by comparing the bottom side of the inner pattern with the top side of the outer pattern
                auto InnerBottomEdge = UniquePatternView::GetPatternEdgeData(PatternDataInner, PatSize, EEdge::Bottom);
                auto OuterTopEdge = UniquePatternView::GetPatternEdgeData(PatternDataOuter, PatSize, EEdge::Top);
                if (InnerBottomEdge == OuterTopEdge)
                {
                    Rules.Top[UniquePatternOuter.Idx].push_back(UniquePatternInner.Idx);
                }

                // Check if the inner pattern can be placed below the outer pattern by comparing the top side of the inner pattern with the bottom side of the outer pattern
                auto InnerTopEdge = UniquePatternView::GetPatternEdgeData(PatternDataInner, PatSize, EEdge::Top);
                auto OuterBottomEdge = UniquePatternView::GetPatternEdgeData(PatternDataOuter, PatSize, EEdge::Bottom);
                if (InnerTopEdge == OuterBottomEdge)
                {
                    Rules.Bottom[UniquePatternOuter.Idx].push_back(UniquePatternInner.Idx);
                }
            }
        }

        return Rules;
    }

    WaveFunctionCollapse::Grid WaveFunctionCollapse::InitializeGrid(const GridSize GridSize, const std::vector<UniquePatternView> &UniquePatterns)
    {
        Grid Grd{};
        Grd.resize(GridSize, std::vector<Cell>(GridSize));

        // Initialize the grid with all patterns as possible in each cell

        // Pre-initialize a "template" cell with all patterns set
        Cell TemplateCell{};
        TemplateCell.PossiblePatterns.reserve(UniquePatterns.size());
        for (PatternIdx i = 0; i < UniquePatterns.size(); ++i)
        {
            TemplateCell.PossiblePatterns.push_back(UniquePatterns[i].Idx);
        }

        // Use the template cell to initialize the grid
        for (auto &Column : Grd)
        {
            std::fill(Column.begin(), Column.end(), TemplateCell);
        }

        return Grd;
    }

    bool WaveFunctionCollapse::Propagate(Grid &Grid, const AdjacencyRules &AdjacencyRules, std::size_t CellX, std::size_t CellY)
    {
        // std::cout << "Propagating constraints...\n";

        // Queue to keep track of cells that need to be processed
        std::queue<std::pair<std::size_t, std::size_t>> Queue;
        Queue.emplace(CellX, CellY);

        // While there are cells in the queue
        while (!Queue.empty())
        {
            auto [CurrentX, CurrentY] = Queue.front();
            Queue.pop();

            Cell &CurrentCell = Grid[CurrentX][CurrentY];

            // Get the possible patterns for the current cell
            std::vector<PatternIdx> CurrentPossiblePatterns;
            if (CurrentCell.CollapsedPatternIdx.has_value())
            {
                // If the cell is collapsed, possible patterns is the single collapsed pattern
                CurrentPossiblePatterns = {CurrentCell.CollapsedPatternIdx.value()};
            }
            else
            {
                // If the cell is not collapsed, use its possible patterns
                CurrentPossiblePatterns = CurrentCell.PossiblePatterns;
            }

            // Define neighbor positions: Left, Right, Top, Bottom
            std::vector<std::pair<std::size_t, std::size_t>> Neighbors;

            // Left neighbor
            if (CurrentX > 0)
                Neighbors.emplace_back(CurrentX - 1, CurrentY);

            // Right neighbor
            if (CurrentX + 1 < Grid.size())
                Neighbors.emplace_back(CurrentX + 1, CurrentY);

            // Top neighbor
            if (CurrentY > 0)
                Neighbors.emplace_back(CurrentX, CurrentY - 1);

            // Bottom neighbor
            if (CurrentY + 1 < Grid[0].size())
                Neighbors.emplace_back(CurrentX, CurrentY + 1);

            // Iterate over each neighbor
            for (auto [NeighborX, NeighborY] : Neighbors)
            {
                Cell &NeighborCell = Grid[NeighborX][NeighborY];

                if (NeighborCell.CollapsedPatternIdx.has_value())
                {
                    // Skip if the neighbor cell is already collapsed
                    continue;
                }

                // Flag to check if the neighbor's possible patterns have changed
                bool Changed = false;

                // New set of possible patterns for the neighbor
                std::vector<PatternIdx> NewPossiblePatterns;

                // Iterate over the neighbor's possible patterns
                for (PatternIdx NeighborPattern : NeighborCell.PossiblePatterns)
                {
                    bool Compatible = false;

                    // Check compatibility with current cell's possible patterns
                    for (PatternIdx Pattern : CurrentPossiblePatterns)
                    {
                        bool IsCompatible = false;

                        // Determine the direction from the current cell to the neighbor
                        if (NeighborX == CurrentX - 1 && NeighborY == CurrentY)
                        {
                            // Neighbor is to the left
                            auto It = AdjacencyRules.Left.find(Pattern);
                            if (It != AdjacencyRules.Left.end())
                            {
                                const auto &AllowedPatterns = It->second;
                                if (std::find(AllowedPatterns.begin(), AllowedPatterns.end(), NeighborPattern) != AllowedPatterns.end())
                                {
                                    IsCompatible = true;
                                }
                            }
                        }
                        else if (NeighborX == CurrentX + 1 && NeighborY == CurrentY)
                        {
                            // Neighbor is to the right
                            auto It = AdjacencyRules.Right.find(Pattern);
                            if (It != AdjacencyRules.Right.end())
                            {
                                const auto &AllowedPatterns = It->second;
                                if (std::find(AllowedPatterns.begin(), AllowedPatterns.end(), NeighborPattern) != AllowedPatterns.end())
                                {
                                    IsCompatible = true;
                                }
                            }
                        }
                        else if (NeighborX == CurrentX && NeighborY == CurrentY - 1)
                        {
                            // Neighbor is above
                            auto It = AdjacencyRules.Top.find(Pattern);
                            if (It != AdjacencyRules.Top.end())
                            {
                                const auto &AllowedPatterns = It->second;
                                if (std::find(AllowedPatterns.begin(), AllowedPatterns.end(), NeighborPattern) != AllowedPatterns.end())
                                {
                                    IsCompatible = true;
                                }
                            }
                        }
                        else if (NeighborX == CurrentX && NeighborY == CurrentY + 1)
                        {
                            // Neighbor is below
                            auto It = AdjacencyRules.Bottom.find(Pattern);
                            if (It != AdjacencyRules.Bottom.end())
                            {
                                const auto &AllowedPatterns = It->second;
                                if (std::find(AllowedPatterns.begin(), AllowedPatterns.end(), NeighborPattern) != AllowedPatterns.end())
                                {
                                    IsCompatible = true;
                                }
                            }
                        }

                        if (IsCompatible)
                        {
                            Compatible = true;
                            break; // No need to check other patterns
                        }
                    }

                    if (Compatible)
                    {
                        // If compatible, keep the neighbor pattern
                        NewPossiblePatterns.push_back(NeighborPattern);
                    }
                    else
                    {
                        // If not compatible, mark as changed
                        Changed = true;
                    }
                }

                if (Changed)
                {
                    NeighborCell.PossiblePatterns = NewPossiblePatterns;

                    if (NeighborCell.PossiblePatterns.empty())
                    {
                        // No possible patterns left: contradiction
                        std::cerr << "Contradiction detected. No possible patterns left for cell (" << NeighborX << ", " << NeighborY << ")\n";
                        return false;
                    }

                    // Add the neighbor cell to the queue to propagate further
                    Queue.emplace(NeighborX, NeighborY);
                }
            }
        }

        return true;
    }

    // Custom hash function for std::pair<std::size_t, std::size_t>
    struct pair_hash
    {
        std::size_t operator()(const std::pair<std::size_t, std::size_t> &p) const
        {
            // A simple combination of the two size_t values.
            // You can use a more sophisticated method if needed.
            return std::hash<std::size_t>()(p.first) ^ (std::hash<std::size_t>()(p.second) << 1);
        }
    };
    bool WaveFunctionCollapse::PropagateWithBackTrack(Grid &Grid, const AdjacencyRules &AdjacencyRules, std::size_t CellX, std::size_t CellY)
    {
        std::cout << "Propagating constraints with backtracking...\n";

        // Stack to keep track of grid states for backtracking
        std::stack<WaveFunctionCollapse::Grid> GridStack;

        // Queue to keep track of cells that need to be processed
        std::queue<std::pair<std::size_t, std::size_t>> Queue;

        // Push the initial state of the grid onto the stack
        GridStack.push(Grid);

        // Start propagation from the given cell
        Queue.emplace(CellX, CellY);

        // Initialize the map to track which patterns have been tried for each cell
        std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t, pair_hash> TriedPatterns;

        while (!Queue.empty())
        {
            auto [CurrentX, CurrentY] = Queue.front();
            Queue.pop();

            Cell &CurrentCell = Grid[CurrentX][CurrentY];

            // Get the possible patterns for the current cell
            std::vector<PatternIdx> CurrentPossiblePatterns;
            if (CurrentCell.CollapsedPatternIdx.has_value())
            {
                // If the cell is collapsed, possible patterns is the single collapsed pattern
                CurrentPossiblePatterns = {CurrentCell.CollapsedPatternIdx.value()};
            }
            else
            {
                // If the cell is not collapsed, use its possible patterns
                CurrentPossiblePatterns = CurrentCell.PossiblePatterns;
            }

            // Define neighbor positions: Left, Right, Top, Bottom
            std::vector<std::pair<std::size_t, std::size_t>> Neighbors;

            // Left neighbor
            if (CurrentX > 0)
                Neighbors.emplace_back(CurrentX - 1, CurrentY);

            // Right neighbor
            if (CurrentX + 1 < Grid.size())
                Neighbors.emplace_back(CurrentX + 1, CurrentY);

            // Top neighbor
            if (CurrentY > 0)
                Neighbors.emplace_back(CurrentX, CurrentY - 1);

            // Bottom neighbor
            if (CurrentY + 1 < Grid[0].size())
                Neighbors.emplace_back(CurrentX, CurrentY + 1);

            // Iterate over each neighbor
            for (auto [NeighborX, NeighborY] : Neighbors)
            {
                Cell &NeighborCell = Grid[NeighborX][NeighborY];

                if (NeighborCell.CollapsedPatternIdx.has_value())
                {
                    // Skip if the neighbor cell is already collapsed
                    continue;
                }

                // New set of possible patterns for the neighbor
                std::vector<PatternIdx> NewPossiblePatterns;

                // Flag to check if the neighbor's possible patterns have changed
                bool Changed = false;

                // Iterate over the neighbor's possible patterns
                for (PatternIdx NeighborPattern : NeighborCell.PossiblePatterns)
                {
                    bool Compatible = false;

                    // Check compatibility with current cell's possible patterns
                    for (PatternIdx Pattern : CurrentPossiblePatterns)
                    {
                        bool IsCompatible = false;

                        // Determine the direction from the current cell to the neighbor
                        if (NeighborX == CurrentX - 1 && NeighborY == CurrentY)
                        {
                            // Neighbor is to the left
                            auto It = AdjacencyRules.Left.find(Pattern);
                            if (It != AdjacencyRules.Left.end())
                            {
                                const auto &AllowedPatterns = It->second;
                                if (std::find(AllowedPatterns.begin(), AllowedPatterns.end(), NeighborPattern) != AllowedPatterns.end())
                                {
                                    IsCompatible = true;
                                }
                            }
                        }
                        else if (NeighborX == CurrentX + 1 && NeighborY == CurrentY)
                        {
                            // Neighbor is to the right
                            auto It = AdjacencyRules.Right.find(Pattern);
                            if (It != AdjacencyRules.Right.end())
                            {
                                const auto &AllowedPatterns = It->second;
                                if (std::find(AllowedPatterns.begin(), AllowedPatterns.end(), NeighborPattern) != AllowedPatterns.end())
                                {
                                    IsCompatible = true;
                                }
                            }
                        }
                        else if (NeighborX == CurrentX && NeighborY == CurrentY - 1)
                        {
                            // Neighbor is above
                            auto It = AdjacencyRules.Top.find(Pattern);
                            if (It != AdjacencyRules.Top.end())
                            {
                                const auto &AllowedPatterns = It->second;
                                if (std::find(AllowedPatterns.begin(), AllowedPatterns.end(), NeighborPattern) != AllowedPatterns.end())
                                {
                                    IsCompatible = true;
                                }
                            }
                        }
                        else if (NeighborX == CurrentX && NeighborY == CurrentY + 1)
                        {
                            // Neighbor is below
                            auto It = AdjacencyRules.Bottom.find(Pattern);
                            if (It != AdjacencyRules.Bottom.end())
                            {
                                const auto &AllowedPatterns = It->second;
                                if (std::find(AllowedPatterns.begin(), AllowedPatterns.end(), NeighborPattern) != AllowedPatterns.end())
                                {
                                    IsCompatible = true;
                                }
                            }
                        }

                        if (IsCompatible)
                        {
                            Compatible = true;
                            break; // No need to check other patterns
                        }
                    }

                    if (Compatible)
                    {
                        // If compatible, keep the neighbor pattern
                        NewPossiblePatterns.push_back(NeighborPattern);
                    }
                    else
                    {
                        // If not compatible, mark as changed
                        Changed = true;
                    }
                }

                if (Changed)
                {
                    // Before making changes, save the current grid state
                    GridStack.push(Grid);

                    NeighborCell.PossiblePatterns = NewPossiblePatterns;

                    if (NeighborCell.PossiblePatterns.empty())
                    {
                        // Contradiction detected
                        std::cerr << "Contradiction detected at cell (" << NeighborX << ", " << NeighborY << "). Backtracking...\n";

                        // Backtrack: Restore the grid from the last saved state
                        if (!GridStack.empty())
                        {
                            Grid = GridStack.top();
                            GridStack.pop();

                            // Identify the cell that caused the contradiction
                            auto cellKey = std::make_pair(NeighborX, NeighborY);
                            auto it = TriedPatterns.find(cellKey);
                            if (it == TriedPatterns.end())
                            {
                                TriedPatterns[cellKey] = 0;
                            }

                            // Get the index of the last tried pattern
                            size_t lastTriedIndex = TriedPatterns[cellKey];

                            if (lastTriedIndex < Grid[NeighborX][NeighborY].PossiblePatterns.size())
                            {
                                // Remove the last tried pattern
                                PatternIdx triedPattern = Grid[NeighborX][NeighborY].PossiblePatterns[lastTriedIndex];
                                Grid[NeighborX][NeighborY].PossiblePatterns.erase(Grid[NeighborX][NeighborY].PossiblePatterns.begin() + lastTriedIndex);
                                std::cerr << "Removing tried pattern " << triedPattern << " from cell (" << NeighborX << ", " << NeighborY << ")\n";
                            }

                            // Increment the tried pattern index
                            TriedPatterns[cellKey]++;

                            // Check if all patterns have been tried
                            if (TriedPatterns[cellKey] >= Grid[NeighborX][NeighborY].PossiblePatterns.size())
                            {
                                // No patterns left to try, need to backtrack further
                                std::cerr << "No patterns left to try for cell (" << NeighborX << ", " << NeighborY << "). Continuing backtrack...\n";
                                Queue.emplace(CellX, CellY);
                            }
                            else
                            {
                                // Assign the next pattern
                                PatternIdx newPattern = Grid[NeighborX][NeighborY].PossiblePatterns[lastTriedIndex];
                                Grid[NeighborX][NeighborY].CollapsedPatternIdx = newPattern;
                                std::cerr << "Assigning new pattern " << newPattern << " to cell (" << NeighborX << ", " << NeighborY << ")\n";
                                // Enqueue the cell to propagate further
                                Queue.emplace(NeighborX, NeighborY);
                            }
                        }
                        else
                        {
                            // No more states to backtrack to
                            std::cerr << "No more states to backtrack to. Failure.\n";
                            return false;
                        }
                    }
                    else
                    {
                        // Add the neighbor cell to the queue to propagate further
                        Queue.emplace(NeighborX, NeighborY);
                    }
                }
            }
        }

        // If the queue is empty and no contradictions occurred, propagation succeeded
        return true;
    }

    std::tuple<std::size_t, std::size_t> WaveFunctionCollapse::ChooseNextCell(const Grid &Grid)
    {
        // Choose the next cell to collapse
        // Typically, the cell with the lowest entropy is chosen

        // Find the cell with the lowest entropy
        std::size_t MinEntropy = std::numeric_limits<std::size_t>::max();
        std::tuple<std::size_t, std::size_t> MinEntropyCell{};

        for (std::size_t y = 0; y < Grid.size(); ++y)
        {
            for (std::size_t x = 0; x < Grid.size(); ++x)
            {
                const auto &Cell = Grid[x][y];

                if (Cell.CollapsedPatternIdx.has_value())
                {
                    // Skip collapsed cells
                    continue;
                }

                const auto Entropy = Cell.PossiblePatterns.size();
                if (Entropy < MinEntropy)
                {
                    MinEntropy = Entropy;
                    MinEntropyCell = {x, y};
                }
            }
        }

        return MinEntropyCell;
    }

    void WaveFunctionCollapse::VisualizeGridEntropy(const Grid &Grid)
    {
        std::cout << "\n\nGrid entropy:\n";

        const auto GridSize = Grid.size();

        // Visualize the grid entropy values
        for (std::size_t y = 0; y < GridSize; ++y)
        {
            for (std::size_t x = 0; x < GridSize; ++x)
            {
                const Cell &Cell = Grid[x][y];
                if (Cell.CollapsedPatternIdx.has_value())
                {
                    std::cout << "X-" << Cell.CollapsedPatternIdx.value() << " ";
                }
                else
                {
                    std::cout << Cell.PossiblePatterns.size() << " ";
                }
            }
            std::cout << "\n";
        }

        std::cout << std::endl;
    }

    void WaveFunctionCollapse::VisualizeAdjacencyRules(const AdjacencyRules &Rules)
    {
        std::cout << "\n\nLeft:\n";
        for (const auto &[Key, Value] : Rules.Left)
        {
            std::cout << Key << ": ";
            for (const auto &V : Value)
            {
                std::cout << V << " ";
            }
            std::cout << "\n";
        }

        std::cout << "\nRight:\n";
        for (const auto &[Key, Value] : Rules.Right)
        {
            std::cout << Key << ": ";
            for (const auto &V : Value)
            {
                std::cout << V << " ";
            }
            std::cout << "\n";
        }

        std::cout << "\nTop:\n";
        for (const auto &[Key, Value] : Rules.Top)
        {
            std::cout << Key << ": ";
            for (const auto &V : Value)
            {
                std::cout << V << " ";
            }
            std::cout << "\n";
        }

        std::cout << "\nBottom:\n";
        for (const auto &[Key, Value] : Rules.Bottom)
        {
            std::cout << Key << ": ";
            for (const auto &V : Value)
            {
                std::cout << V << " ";
            }
            std::cout << "\n";
        }

        std::cout << std::endl;
    }

    bool WaveFunctionCollapse::Run(const RunParameters &Params)
    {
        // Load the input image
        auto InputImage = LoadImage(Params.InputImageFilename);

        // Extract unique patterns from the input image
        auto UniquePatterns = ExtractPatterns(InputImage, Params.PatSize);

        if (UniquePatterns.empty())
        {
            return false;
        }

        // Extract adjacency rules from the patterns
        const auto AdjacencyRules = ExtractAdjacencyRules(UniquePatterns, Params.PatSize, InputImage);

        // Visualize the adjacency rules
        // VisualizeAdjacencyRules(AdjacencyRules);

        // Initialize the grid with all patterns as possible in each cell
        auto Grid = InitializeGrid(Params.GridSize, UniquePatterns);

        // Visualize the initial grid entropy values
        // VisualizeGridEntropy(Grid);

        // Choose a random cell to collapse
        std::random_device RD{};
        std::mt19937 RNG{RD()};
        std::uniform_int_distribution<std::size_t> Dist{0, Params.GridSize - 1};

        const auto MaxRetries{1000};
        std::size_t RetryCount{0};

        while (true)
        {
            auto [RandomX, RandomY] = ChooseNextCell(Grid);
            auto &CellToCollapse = Grid[RandomX][RandomY];

            // Collapse the random cell
            CellToCollapse.CollapsedPatternIdx = CellToCollapse.PossiblePatterns[Dist(RNG) % CellToCollapse.PossiblePatterns.size()];

            // Propagate constraints
            if (!Propagate(Grid, AdjacencyRules, RandomX, RandomY))
            {
                if (++RetryCount >= MaxRetries)
                {
                    std::cerr << "Max retries reached. Aborting...\n";

                    VisualizeGridEntropy(Grid);

                    return false;
                }

                // Reset the grid
                Grid = InitializeGrid(Params.GridSize, UniquePatterns);
                continue;
            }

            if (std::all_of(Grid.begin(), Grid.end(), [](const auto &Row)
                            { return std::all_of(Row.begin(), Row.end(), [](const auto &Cell)
                                                 { return Cell.CollapsedPatternIdx.has_value(); }); }))
            {
                // All cells are collapsed
                break;
            }
        }

        // Visualize the grid entropy values
        // VisualizeGridEntropy(Grid);

        // Save the output image
        const auto NumPatternsPerRow = ((std::get<0>(InputImage) - Params.PatSize) / (Params.PatSize - 1)) + 1;
        const auto NumPatternsPerCol = ((std::get<1>(InputImage) - Params.PatSize) / (Params.PatSize - 1)) + 1;
        const auto OutputWidth = Params.GridSize * Params.PatSize;
        const auto OutputHeight = Params.GridSize * Params.PatSize;

        std::vector<std::uint8_t> OutputImage(OutputWidth * OutputHeight * 4);

        for (std::size_t y = 0; y < Params.GridSize; ++y)
        {
            for (std::size_t x = 0; x < Params.GridSize; ++x)
            {
                const auto &Cell = Grid[x][y];

                const auto [PatternX, PatternY] = UniquePatternView::ComputePatternCoord(Cell.CollapsedPatternIdx.value(), NumPatternsPerRow, NumPatternsPerCol);
                const auto PatternData = UniquePatternView::GetPatternData({PatternX, PatternY}, Params.PatSize, InputImage);

                for (std::size_t py = 0; py < Params.PatSize; ++py)
                {
                    for (std::size_t px = 0; px < Params.PatSize; ++px)
                    {
                        // Get the pixel data from the pattern
                        auto &Pixel = PatternData[py * Params.PatSize + px];

                        // std::cout << "Pixel: " << (int)std::get<0>(Pixel) << " " << (int)std::get<1>(Pixel) << " " << (int)std::get<2>(Pixel) << " " << (int)std::get<3>(Pixel) << "\n";

                        // Calculate the index in the output image
                        std::size_t OutputIndex = ((y * Params.PatSize + py) * OutputWidth + (x * Params.PatSize + px)) * 4;

                        // Copy the pixel data to the output image
                        OutputImage[OutputIndex] = std::get<0>(Pixel);
                        OutputImage[OutputIndex + 1] = std::get<1>(Pixel);
                        OutputImage[OutputIndex + 2] = std::get<2>(Pixel);
                        OutputImage[OutputIndex + 3] = std::get<3>(Pixel);
                    }
                }
            }
        }

        stbi_write_png(Params.OutputImageFilename.c_str(), OutputWidth, OutputHeight, 4, OutputImage.data(), OutputWidth * 4);

        return true;
    }
}