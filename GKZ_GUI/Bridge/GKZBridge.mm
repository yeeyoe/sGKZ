#import "GKZBridge.h"

#import "../CoreAdapter/GuiSolverAdapter.hpp"
#include <atomic>
#include <memory>
#import <objc/runtime.h>

@interface GKZGuiResult ()
@property(nonatomic) BOOL certified;
@property(nonatomic) BOOL cancelled;
@property(nonatomic, copy) NSString *message;
@property(nonatomic, copy) NSArray<NSString *> *sigmaExact;
@property(nonatomic, copy) NSArray<NSString *> *sigmaVeeExact;
@property(nonatomic, copy) NSArray<NSArray<NSNumber *> *> *triangulationFaces;
@property(nonatomic, copy) NSArray<NSArray<NSValue *> *> *subdivisionCells;
@property(nonatomic) NSInteger iterations;
@property(nonatomic) NSUInteger activeSize;
@end

@interface GKZCancellationToken ()
- (std::atomic_bool *)atomicFlag;
@end

@implementation GKZCancellationToken
static char GKZCancellationFlagKey;
- (instancetype)init {
  self = [super init];
  if (self) {
    auto *flag = new std::atomic_bool(false);
    objc_setAssociatedObject(self, &GKZCancellationFlagKey,
                             [NSValue valueWithPointer:flag],
                             OBJC_ASSOCIATION_RETAIN_NONATOMIC);
  }
  return self;
}
- (void)dealloc { delete [self atomicFlag]; }
- (void)cancel { [self atomicFlag]->store(true); }
- (BOOL)isCancelled { return [self atomicFlag]->load(); }
- (std::atomic_bool *)atomicFlag {
  NSValue *value = objc_getAssociatedObject(self, &GKZCancellationFlagKey);
  return static_cast<std::atomic_bool *>(value.pointerValue);
}
@end

@implementation GKZGuiResult
@end

@implementation GKZBridge
+ (GKZGuiResult *)solvePoints:(NSArray<NSValue *> *)points
                     progress:(GKZProgressBlock)progress
                        token:(GKZCancellationToken *)token {
  std::vector<gkz_gui::GuiPoint> input;
  input.reserve(points.count);
  for (NSValue *value in points) {
    CGPoint point = value.pointValue;
    input.push_back({static_cast<long long>(point.x), static_cast<long long>(point.y)});
  }
  const auto result = gkz_gui::solve_points(
      input, [token atomicFlag],
      [progress, token](int iterations, std::size_t active, const std::string& status) {
        if (token.cancelled) return;
        NSString *statusString = [NSString stringWithUTF8String:status.c_str()] ?: @"";
        if (progress) dispatch_async(dispatch_get_main_queue(), ^{
          progress(iterations, active, statusString);
        });
      });
  auto output = [GKZGuiResult new];
  output.certified = result.certified;
  output.cancelled = result.cancelled || token.cancelled;
  output.message = [NSString stringWithUTF8String:result.message.c_str()] ?: @"";
  NSMutableArray *sigma = [NSMutableArray array];
  for (const auto& item : result.sigma_exact)
    [sigma addObject:[NSString stringWithUTF8String:item.c_str()] ?: @""];
  output.sigmaExact = sigma;
  NSMutableArray *sigmaVee = [NSMutableArray array];
  for (const auto& item : result.sigma_vee_exact)
    [sigmaVee addObject:[NSString stringWithUTF8String:item.c_str()] ?: @""];
  output.sigmaVeeExact = sigmaVee;
  NSMutableArray *faces = [NSMutableArray array];
  for (const auto& face : result.triangulation_faces)
    [faces addObject:@[@(face[0]), @(face[1]), @(face[2])]];
  output.triangulationFaces = faces;
  NSMutableArray *cells = [NSMutableArray array];
  for (const auto& cell : result.subdivision_cells) {
    NSMutableArray *vertices = [NSMutableArray array];
    for (const auto& point : cell)
      [vertices addObject:[NSValue valueWithPoint:NSMakePoint(point.first, point.second)]];
    [cells addObject:vertices];
  }
  output.subdivisionCells = cells;
  output.iterations = result.iterations;
  output.activeSize = result.active_size;
  return output;
}
@end
