abstract type AbstractFunction end

struct SimpleFunction <: AbstractFunction 
    f::Function
end

struct PartionedFunction <: AbstractFunction
    f1::Function
    f2::Function
end
