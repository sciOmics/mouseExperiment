// Mobile device incompatibility overlay
// Only injects the overlay when the viewport is <= 768 px.
$(document).ready(function() {
  if (window.innerWidth > 768) return;

  var html =
    '<div id="mobile-overlay" style="' +
      'position:fixed;top:0;left:0;right:0;bottom:0;z-index:99999;' +
      'background:#f5f5f5;color:#333;display:flex;flex-direction:column;' +
      'align-items:center;justify-content:center;text-align:center;' +
      'padding:2rem 1.5rem;overflow-y:auto;">' +
      '<div style="font-size:4rem;margin-bottom:1rem;">&#x1F4BB;</div>' +
      '<h2 style="font-size:1.4rem;margin-bottom:0.75rem;">Desktop Recommended</h2>' +
      '<p style="font-size:1rem;line-height:1.5;max-width:360px;' +
        'margin:0 auto 1.5rem;color:#555;">This dashboard requires a larger screen ' +
        'for the best experience. Please visit on a desktop or laptop computer.</p>' +
      '<button id="mobile-dismiss-btn" style="' +
        'padding:0.5rem 1.5rem;border:1px solid #999;border-radius:4px;' +
        'background:transparent;color:#555;font-size:0.85rem;cursor:pointer;' +
      '">Continue anyway</button>' +
    '</div>';

  var $overlay = $(html);
  $('body').append($overlay);

  $overlay.on('click', '#mobile-dismiss-btn', function() {
    $overlay.remove();
  });
});
