module ZeissMicroscopyFormat

using Mmap
using Dates

using FileIO
using OMETIFF
using EzXML
using FixedPointNumbers
using ColorTypes
using ImageMetadata
using Unitful
using LazyStack

# Documentation on this format is proprietary but available by request for free:
#    https://www.zeiss.com/microscopy/int/products/microscope-software/zen/czi.html
# The license bars posting the API documentation online, but it does not prohibit other usages.

# Segment IDs
segment_id(str) = FileIO.padzeros(str, 16)
const RAWDIR   = segment_id("ZISRAWDIRECTORY")
const SUBBLOCK = segment_id("ZISRAWSUBBLOCK")
const METADATA = segment_id("ZISRAWMETADATA")
const ATTACH   = segment_id("ZISRAWATTACH")
const ATTDIR   = segment_id("ZISRAWATTDIR")
const DELETED  = segment_id("DELETED")

const pixeltypes = Dict{Int32,DataType}(
    0 => Gray{N0f8},
    1 => Gray{N0f16},
    2 => Gray{Float32},
    3 => BGR{N0f8},
    4 => BGR{N0f16},
    8 => BGR{Float32},
    9 => BGRA{N0f8},
    # unsupported: 10 => Gray{Complex{Float32}}, 11 => BGR{Complex{Float32}}
    12 => Gray{N0f32}
)

const zeissdtfmt = dateformat"y-m-dTH:M:S.s"

struct DimensionEntryDV
    dimension::Char
    start::Int32
    size::Int32
    startcoord::Float32
    storedsize::Int32
end
# DimensionEntryDV(s::Stream) = reinterpret(DimensionEntryDV, read(s, sizeof(DimensionEntryDV)))[]
DimensionEntryDV(s::Stream) = DimensionEntryDV(
    Char(read(s, UInt32)),
    read(s, Int32),
    read(s, Int32),
    read(s, Float32),
    read(s, Int32),
)

function Base.show(io::IO, dime::DimensionEntryDV)
    print(io, "Dim")
    if !iszero(dime.storedsize)
        print(io, '*')
    end
    print(io, ' ', dime.dimension, ' ', dime.size, " at ", dime.startcoord)
end

struct DirectoryEntryDVFixed
    pixeltype::Int32
    fileposition::Int64
    filepart::Int32
    compression::Int32
    pyramidtype::UInt8
    spare::UInt8
    spare4::UInt32
    dimensioncount::Int32
end

struct DirectoryEntryDV
    fixed::DirectoryEntryDVFixed
    dimensionentries::Vector{DimensionEntryDV}
end

function DirectoryEntryDV(s::Stream)
    c1, c2 = read(s, UInt8), read(s, UInt8)
    @assert Char(c1) == 'D' && Char(c2) == 'V'
    # fixed = reinterpret(DirectoryEntryDVFixed, read(s, sizeof(DirectoryEntryDVFixed)))[]
    fixed = DirectoryEntryDVFixed(read(s, Int32),
                                  read(s, Int64),
                                  read(s, Int32),
                                  read(s, Int32),
                                  read(s, UInt8),
                                  read(s, UInt8),    # spare
                                  read(s, UInt32),   # spare4
                                  read(s, Int32),
    )
    @assert iszero(fixed.compression)   # FIXME: support compression
    # entries = reinterpret(DimensionEntryDV, read(s, fixed.dimensioncount * sizeof(DimensionEntryDV)))
    entries = [DimensionEntryDV(s) for _ = 1:fixed.dimensioncount]
    return DirectoryEntryDV(fixed, entries)
end
function Base.getproperty(d::DirectoryEntryDV, name::Symbol)
    name == :dimensionentries && return getfield(d, :dimensionentries)
    return getproperty(getfield(d, :fixed), name)
end

Base.size(d::DirectoryEntryDV) = ntuple(d.dimensioncount) do i
    Int(d.dimensionentries[i].size)
end
Base.names(d::DirectoryEntryDV) = ntuple(d.dimensioncount) do i
    Symbol(d.dimensionentries[i].dimension)
end

struct SubBlockSegmentSizes
    metasize::Int32
    attachsize::Int32
    datasize::Int64
end
struct DataSegment{C<:Colorant,N}
    offset::Int64
    dims::NTuple{N,Int}
    dimnames::NTuple{N,Char}
    metadata
end

function read_subblock(s::Stream, d::DirectoryEntryDV) where C<:Colorant
    pos = position(s)
    sbss = reinterpret(SubBlockSegmentSizes, read(s, sizeof(SubBlockSegmentSizes)))[1]
    seek(s, pos+256)
    md = parsexml(String(read(s, sbss.metasize)))
    return DataSegment{pixeltypes[d.pixeltype],Int(d.dimensioncount)}(position(s), )
end

function load(f::File{format"CZI"})
    open(f) do s
        skipmagic(s)  # skip over the magic bytes
        return load(s)
    end
end

function load(s::Stream{format"CZI"}; keywords...)
    # s is already positioned after the magic bytes
    allocated_size = read(s, Int64)
    used_size = read(s, Int64)
    major = read(s, Int32)
    @assert major == 1
    minor = read(s, Int32)
    @assert minor == 0
    read(s, Int32)  # reserved
    read(s, Int32)  # reserved
    primaryfileguid = read(s, UInt128)
    fileguid = read(s, UInt128)
    filepart = read(s, Int32)
    dirpos  = read(s, Int64)
    metapos = read(s, Int64)
    updatepending = read(s, Int32)
    attpos = read(s, Int64)
    @assert iszero(updatepending) "update is pending, please try later"
    @assert iszero(filepart) "multi-part files not yet supported"

    # Metadata
    seek(s, metapos)
    sid = read(s, 16)
    @assert sid == METADATA
    read(s, 16)
    xml_size, attach_size = read(s, Int32), read(s, Int32)
    seek(s, metapos + 256 + 32)
    rawxml = String(read(s, xml_size))
    # idx = findfirst(==('\0'), rawxml)
    # if idx !== nothing
    #     rawxml = rawxml[1:prevind(rawxml, idx)]
    # end
    omexml = OMETIFF.root(OMETIFF.parsexml(rawxml))

    # Subblocks
    seek(s, dirpos)
    sid = read(s, 16)
    @assert sid == RAWDIR
    read(s, 16)
    entry_count = read(s, Int32)
    seek(s, dirpos + 128 + 32)
    entries = [DirectoryEntryDV(s) for i = 1:entry_count]
    pixeltype = first(entries).pixeltype
    nd = first(entries).dimensioncount
    sz = size(first(entries))
    @show names(first(entries))
    # @show sz map(size, entries)
    @assert all(d -> d.pixeltype == pixeltype, entries) "All pixel types must be identical"
    @assert all(d -> d.dimensioncount == nd, entries) "The number of dimensions must be consistent"
    @assert all(d -> size(d)[1:end-1] == sz[1:end-1], entries) "Leading sizes must be identical"

    # Data
    # FIXME? To avoid another layer of `ReinterpretArray`, it's easiest to interpret the raw data
    # as having the eltype of the image. However, it's possible the chunk boundaries will not
    # always be aligned commensurately. For now let's run with this, but it may require generalization.
    T = pixeltypes[pixeltype]
    nb = sizeof(T)::Int
    @assert all(d -> iszero(d.fileposition % nb), entries) "FIXME: chunk boundaries are not aligned"
    C, axs, shear = layout(omexml, T)
    @show C axs shear
    seek(s, 0)
    mm = Mmap.mmap(s.io, Vector{T}, filesize(s.io) ÷ nb)
    blocks = [makeview(mm, entry, nb) for entry in entries]

    return ImageMeta(stack(blocks); xml=omexml, entries)  # FIXME
end

function makeview(v, entry, nb)
    start = entry.fileposition ÷ nb
    sz = size(entry)
    nel = prod(sz)
    # Drop the color channel ('C') if eltype(v) <: Gray
    if eltype(v) <: AbstractGray
        nms = names(entry)
        idx = findfirst(==(:C), nms)
        if idx !== nothing
            sz = (sz[1:idx-1]..., sz[idx+1:end]...)
        end
    end
    return reshape(view(v, start:start+nel-1), sz)
end

function layout(omexml, ::Type{T}) where T
    imagenode = only(filter!(node -> !isempty(eachelement(node)), findall("//Image", omexml)))
    axs = Dict{Symbol,Any}()
    szs = Dict{Symbol,Int}()
    wavelengths = Float32[]
    shear = nothing
    for node in eachelement(imagenode)
        nn = nodename(node)
        if startswith(nn, "Size")
            sym = Symbol(nn[end])
            szs[sym] = parse(Int, nodecontent(node))
        elseif nn == "Dimensions"
            for subnode in eachelement(node)
                nn = nodename(subnode)
                if nn == "Channels"
                    for cnode in eachelement(subnode)
                        for anode in eachelement(cnode)
                            nn = nodename(anode)
                            if nn == "EmissionWavelength"
                                push!(wavelengths, parse(Float32, nodecontent(anode)))
                            end
                        end
                    end
                elseif nn == "T"
                    starttime, start, increment = nothing, nothing, nothing
                    for anode in eachelement(subnode)
                        nn = nodename(anode)
                        if nn == "StartTime"
                            starttime = nodecontent(anode)
                        elseif nn == "Positions"
                            start, increment = parse_positions(anode)
                        end
                    end
                    idx = findfirst('.', starttime)
                    starttime = DateTime(starttime[1:idx+3], zeissdtfmt)   # TODO: do something with this?
                    axs[:T] = range(start * u"s", step=increment * u"s", length=szs[:T])
                elseif nn == "Z"
                    start, increment = nothing, nothing
                    for anode in eachelement(subnode)
                        nn = nodename(anode)
                        if nn == "ZAxisShear"
                            shear = nodecontent(anode)
                        elseif nn == "Positions"
                            start, increment = parse_positions(anode)
                        end
                    end
                    axs[:Z] = range(start * u"μm", step=increment * u"μm", length=szs[:Z])
                end
            end
        end
    end
    return wavelengths, axs, shear
end

function parse_positions(node)
    node = only(eachelement(node))  # Interval
    startnode = firstelement(node)
    @assert nodename(startnode) == "Start"
    start = nodecontent(startnode)
    incnode = nextelement(startnode)
    @assert nodename(incnode) == "Increment"
    increment = nodecontent(incnode)
    return parse(Float64, start), parse(Float64, increment)
end

end
